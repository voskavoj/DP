import itertools

from astropy.coordinates import EarthLocation, ITRS
from astropy.time import Time, TimeDelta
import astropy.units as unit
import numpy as np
import matplotlib.pyplot as plt

from src.navigation.calculations import m_to_deg_lat, m_to_deg_lon, latlon_distance
from src.satellites.predictions import predict_satellite_positions
from src.navigation.data_processing import nav_data_to_array, find_curves
from src.navigation.data_processing import NavDataArrayIndices as IDX
from src.config.parameters import CurveFitMethodParameters
from src.config.locations import LOCATIONS
from src.utils.data import dump_data
from src.utils.plots import plot_results_of_iterative_position_finding, plot_analyzed_curve, \
    plot_measured_vs_trial_curve
from src.config.setup import DEBUG


# Constants
C = 299792458  # m/s
LON_HOME, LAT_HOME, ALT_HOME = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]


class IterationResults:
    first = 0
    found = 1
    limit = 2


def estimate_zero_doppler_shift_position(detected_curves, satellites):
    """
        Estimate the position of the satellite at the moment of zero doppler shift, as a weighted average of all curves

    :param detected_curves: list of curves
    :param satellites: satellites
    :return: lat, lon; or 0, 0 if no valid curves were found
    """

    estimated_init_locations = list()
    for curve_array in detected_curves:
        sat = satellites[str(int(curve_array[0, IDX.sat_id]))]
        if DEBUG.log_detected_curves:
            print(f"\t {sat.name}, {curve_array[0, IDX.fb]:.0f} Hz, {curve_array.shape[0]: 5d} samples")

        # interpolate moment of zero doppler shift
        dopp_shift_arr = curve_array[:, IDX.f] - curve_array[:, IDX.fb]
        # for interpolation the doppler shift must be increasing -> flip if necessary
        if dopp_shift_arr[0] > dopp_shift_arr[-1]:
            dopp_shift_arr = dopp_shift_arr[::-1]
        pass_time = np.interp(0, xp=dopp_shift_arr, fp=curve_array[:, IDX.t], left=0, right=0)

        # if interpolation failed, skip the curve
        if pass_time == 0:
            continue

        # predict satellite position at pass time
        pass_pos = predict_satellite_positions([sat.satrec], Time(np.array([pass_time]), format="unix"))[0, 0]
        pass_pos_ground = pass_pos.earth_location.geodetic

        estimated_init_locations.append((pass_pos_ground.lat.value, pass_pos_ground.lon.value, curve_array.shape[0]))

    # if no valid points were found, return 0, 0 as a last ditch solution
    if not estimated_init_locations:
        return 0, 0

    # calculate lat, lon as a weighted average (weights are the length of the curve)
    sum_lat = sum([lat * weight for lat, _, weight in estimated_init_locations])
    sum_lon = sum([lon * weight for _, lon, weight in estimated_init_locations])
    sum_weight = sum([weight for _, _, weight in estimated_init_locations])

    return sum_lat / sum_weight, sum_lon / sum_weight


def check_trial_curve(lat, lon, alt, off, dft, measured_curve, r_sat_arr, v_sat_arr, plot=False):
    curve_len = measured_curve.shape[0]

    trial_curve = np.empty((curve_len, 2))
    trial_curve[:, 0] = measured_curve[:, 0]  # times are the same
    r_user_arr = (EarthLocation.from_geodetic(lon, lat, alt)
                  .get_itrs().cartesian.without_differentials())

    vs, rs, ru = v_sat_arr.T * 1000, r_sat_arr.T * 1000, r_user_arr * np.ones(curve_len)
    ru = np.array([ru.x.to("km").value, ru.y.to("km").value, ru.z.to("km").value]) * 1000
    f_b = measured_curve[:, 2]

    # Calculate range rate: Ro_dot = (V_s - V_u) * (r_s - r_u) / ||r_s - r_u||
    rel_vel = np.sum(vs * (rs - ru) / np.linalg.norm(rs - ru, axis=0), axis=0)
    # Calculate doppler frequency
    f_d = -1 * rel_vel * f_b / C

    # Calculate drift
    f_drift = (trial_curve[:, 0] - np.min(trial_curve[:, 0])) * dft

    # Construct trial curve
    trial_curve[:, 1] = f_d + off + f_drift

    # Calculate variance
    sum_of_squares = np.sum((measured_curve[:, 1] - trial_curve[:, 1]) ** 2)

    if plot:
        plot_measured_vs_trial_curve(measured_curve, trial_curve, lat, lon, alt, off)

    return sum_of_squares


def _adjust_step(step, diff, step_limit, force=False):
    if (force or diff == 0) and step > step_limit:
        step = max(step_limit, round(step * 0.75))
    return step


class RatioCycleDetector:
    """
    Detects if the algorithm is stuck in a cycle
    """
    def __init__(self):
        self.last_ratios = None
        self.second_to_last_ratios = None
        self.was_cycle = False

    def update(self, ratios):
        self.was_cycle = False
        self.second_to_last_ratios = self.last_ratios
        self.last_ratios = ratios

    def is_cycle(self):
        if self.last_ratios is None or self.second_to_last_ratios is None:
            return False

        cycle = [int(self.last_ratios[i] + self.second_to_last_ratios[i]) for i in range(len(self.last_ratios))] == [0] * len(self.last_ratios)

        if cycle:
            self.last_ratios, self.second_to_last_ratios = None, None

        self.was_cycle = cycle
        return cycle


def fit_curve(results, lat_0, lon_0, alt_0, off_0, dft_0, measured_curve, r_sat_arr, v_sat_arr, params):
    """
    Iterative curve fitting method

    "sos" is the metric

    :param results: list of results
    :param lat_0: initial latitude
    :param lon_0: initial longitude
    :param alt_0: initial altitude
    :param off_0: initial offset
    :param dft_0: initial drift
    :param measured_curve: array of measured nav data
    :param r_sat_arr: array of satellite positions
    :param v_sat_arr: array of satellite velocities
    :param params: parameters
    :return: final metric, lat, lon, alt, off, dft, results
    """
    iteration_result = IterationResults.limit
    cycle_detector = RatioCycleDetector()

    # initialization
    step_lat = params.lat.init_step
    step_lon = params.lon.init_step
    step_alt = params.alt.init_step
    step_off = params.off.init_step
    step_dft = params.dft.init_step
    lat = lat_0
    lon = lon_0
    alt = alt_0
    off = off_0
    dft = dft_0

    # initial SOS
    sos = check_trial_curve(lat, lon, alt, off, 0, measured_curve, r_sat_arr, v_sat_arr)
    results.append([sos, -1, lat, lon, alt, off, dft])

    for i in range(params.iteration_limit):
        # nORTH, sOUTH, wEST, eAST, uP, dOWN, mORE_OFFSET, lESS_OFFSET

        # 1. calculate SOS for each direction
        lat_n = lat + m_to_deg_lat(step_lat * 1, lat)
        sos_n = check_trial_curve(lat_n, lon, alt, off, dft, measured_curve, r_sat_arr, v_sat_arr)

        lat_s = lat + m_to_deg_lat(step_lat * -1, lat)
        sos_s = check_trial_curve(lat_s, lon, alt, off, dft, measured_curve, r_sat_arr, v_sat_arr)

        lon_w = lon + m_to_deg_lon(step_lon * -1, lat)
        sos_w = check_trial_curve(lat, lon_w, alt, off, dft, measured_curve, r_sat_arr, v_sat_arr)

        lon_e = lon + m_to_deg_lon(step_lon * 1, lat)
        sos_e = check_trial_curve(lat, lon_e, alt, off, dft, measured_curve, r_sat_arr, v_sat_arr)

        alt_u = alt + step_alt * 1
        sos_u = check_trial_curve(lat, lon, alt_u, off, dft, measured_curve, r_sat_arr, v_sat_arr)

        alt_d = alt + step_alt * -1
        sos_d = check_trial_curve(lat, lon, alt_d, off, dft, measured_curve, r_sat_arr, v_sat_arr)

        off_m = off + step_off * 1
        sos_m = check_trial_curve(lat, lon, alt, off_m, dft, measured_curve, r_sat_arr, v_sat_arr)

        off_l = off + step_off * -1
        sos_l = check_trial_curve(lat, lon, alt, off_l, dft, measured_curve, r_sat_arr, v_sat_arr)

        dft_m = dft + step_dft * 1
        sos_dm = check_trial_curve(lat, lon, alt, off, dft_m, measured_curve, r_sat_arr, v_sat_arr)

        dft_l = dft + step_dft * -1
        sos_dl = check_trial_curve(lat, lon, alt, off, dft_l, measured_curve, r_sat_arr, v_sat_arr)

        # 2. calculate differences
        sos_n = max(0, sos - sos_n)
        sos_s = max(0, sos - sos_s)
        sos_w = max(0, sos - sos_w)
        sos_e = max(0, sos - sos_e)
        sos_u = max(0, sos - sos_u)
        sos_d = max(0, sos - sos_d)
        sos_m = max(0, sos - sos_m)
        sos_l = max(0, sos - sos_l)
        sos_dm = max(0, sos - sos_dm)
        sos_dl = max(0, sos - sos_dl)

        diff_lat = sos_n - sos_s
        diff_lon = sos_e - sos_w
        diff_alt = sos_u - sos_d
        diff_off = sos_m - sos_l
        diff_dft = sos_dm - sos_dl

        # 3. check if current position is the best

        if (alt <= params.alt.lower_bound and diff_alt < 0) or (alt >= params.alt.upper_bound and diff_alt > 0):
            diff_alt = 0

        step_lat = _adjust_step(step_lat, diff_lat, params.lat.step_limit)
        step_lon = _adjust_step(step_lon, diff_lon, params.lon.step_limit)
        step_alt = _adjust_step(step_alt, diff_alt, params.alt.step_limit)
        step_off = _adjust_step(step_off, diff_off, params.off.step_limit)
        step_dft = _adjust_step(step_dft, diff_dft, params.dft.step_limit)

        if (diff_lat == 0 and diff_lon == 0 and diff_off == 0 and diff_alt == 0 and diff_dft == 0
                and step_lat <= params.lat.step_limit and step_lon <= params.lon.step_limit):

            results.append([sos, i, lat, lon, alt, off, dft])
            print(f"Iteration {i + 1:03d}: lat {lat:03.3f}°, lon {lon:02.3f}°, alt {alt:04.0f} m, off {off:05.0f} Hz, dft {dft:05.3f} Hz/s, "
                  f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m, "
                  f"ITERATION END")
            if i == 0:
                iteration_result = IterationResults.first
            else:
                iteration_result = IterationResults.found
            break

        # 4. calculate ratios
        ratio_lat = diff_lat / abs(diff_lat) if diff_lat else 0
        ratio_lon = diff_lon / abs(diff_lon) if diff_lon else 0
        ratio_alt = diff_alt / abs(diff_alt) if diff_alt else 0
        ratio_off = diff_off / abs(diff_off) if diff_off else 0
        ratio_dft = diff_dft / abs(diff_dft) if diff_dft else 0

        # 5. move trial position and calculate new SOS
        # lat = lat + m_to_deg_lat(step_lat * ratio_lat, lat)
        # lon = lon + m_to_deg_lon(step_lon * ratio_lon, lat)
        # alt = min(params.alt.upper_bound, max(params.alt.lower_bound, alt + step_alt * ratio_alt))
        # off = off + step_off * ratio_off
        # dft = dft + step_dft * ratio_dft

        # check if the algorithm is not stuck in a cycle
        cycle_detector.update([ratio_lat, ratio_lon, ratio_alt, ratio_off])
        if cycle_detector.is_cycle():  # cycle detection
            step_lat = _adjust_step(step_lat, diff_lat, params.lat.step_limit, force=ratio_lat != 0)
            step_lon = _adjust_step(step_lon, diff_lon, params.lon.step_limit, force=ratio_lon != 0)
            step_alt = _adjust_step(step_alt, diff_alt, params.alt.step_limit, force=ratio_alt != 0)
            step_off = _adjust_step(step_off, diff_off, params.off.step_limit, force=ratio_off != 0)
            step_dft = _adjust_step(step_dft, diff_dft, params.dft.step_limit, force=ratio_dft != 0)

        lat = lat + m_to_deg_lat(step_lat * ratio_lat, lat)
        lon = lon + m_to_deg_lon(step_lon * ratio_lon, lat)
        alt = min(params.alt.upper_bound, max(params.alt.lower_bound, alt + step_alt * ratio_alt))
        off = off + step_off * ratio_off
        dft = dft + step_dft * ratio_dft

        sos = check_trial_curve(lat, lon, alt, off, dft, measured_curve, r_sat_arr, v_sat_arr)

        # 6. log results
        results.append([sos, i, lat, lon, alt, off, dft])
        if DEBUG.log_detail_progress:
            print(f"\tIteration {i + 1:03d}: lat {lat:05.3f}°, lon {lon:05.3f}°, alt {alt:04.0f} m, off {off:05.0f} Hz, dft {dft:05.3f} Hz/s, SOS {sos:09.0f}, "
                  f"step_lat {step_lat:05.0f} m, step_lon {step_lon:05.0f} m, step_alt {step_alt:03.0f} m, step_off {step_off:04.0f} Hz, step_dft {step_dft:05.3f} Hz/s,"
                  f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m, "
                  f"rlat {ratio_lat: 2.0f}, rlon {ratio_lon: 2.0f}, ralt {ratio_alt: 2.0f}, roff {ratio_off: 2.0f}, rdft {ratio_dft: 2.0f}{', Cycle' if cycle_detector.was_cycle else ''}")

    return iteration_result, lat, lon, alt, off, dft, results


def _generate_grid_search_steps(val, step, params):
    """
    Generate grid search steps for value (state vector member)
    :param val: state vector member
    :param step: step
    :param params: parameters
    :return: steps
    """
    if step == 0:
        return [val]

    steps = list()

    for i in range(-int(params.grid_size/2), int(params.grid_size/2) + 1):
        new_val = val + step * i
        if params.lower_bound is not None and new_val < params.lower_bound:
            continue
        if params.upper_bound is not None and new_val > params.upper_bound:
            continue

        steps.append(new_val)
    if not steps:
        steps = [val]
    return steps


def fit_curve_grid(results, lat_0, lon_0, alt_0, off_0, dft_0, measured_curve, r_sat_arr, v_sat_arr, params):
    """
    Grid-based curve fitting method

    "sos" is the metric

    :param results: list of results
    :param lat_0: initial latitude
    :param lon_0: initial longitude
    :param alt_0: initial altitude
    :param off_0: initial offset
    :param dft_0: initial drift
    :param measured_curve: array of measured nav data
    :param r_sat_arr: array of satellite positions
    :param v_sat_arr: array of satellite velocities
    :param params: parameters
    :return: result, final lat, lon, alt, off, dft, results
    """

    lat = lat_0
    lon = lon_0
    alt = alt_0
    off = off_0
    dft = dft_0

    res = list()

    for i in range(params.grid_search_depth):
        step_lat = m_to_deg_lat(params.lat.steps[i], lat)
        step_lon = m_to_deg_lat(params.lon.steps[i], lat)
        step_alt = params.alt.steps[i]
        step_off = params.off.steps[i]
        step_dft = params.dft.steps[i]

        lats = _generate_grid_search_steps(lat, step_lat, params.lat)
        lons = _generate_grid_search_steps(lon, step_lon, params.lon)
        alts = _generate_grid_search_steps(alt, step_alt, params.alt)
        offs = _generate_grid_search_steps(off, step_off, params.off)
        dfts = _generate_grid_search_steps(dft, step_dft, params.dft)

        for lat, lon, alt, off, dft in itertools.product(lats, lons, alts, offs, dfts):
            sos = check_trial_curve(lat, lon, alt, off, dft, measured_curve, r_sat_arr, v_sat_arr)
            res.append([sos, 0, lat, lon, alt, off, dft])

        sos, _, lat, lon, alt, off, dft = min(res, key=lambda x: x[0])
        res.append([sos, 0, lat, lon, alt, off, dft])
        if DEBUG.log_detail_progress:
            print(f"Grid {i + 1:03d}: lat {lat:05.3f}°, lon {lon:05.3f}°, alt {alt:04.0f} m, off {off:05.0f} SOS {sos:09.0f}, "
                  f"step_lat {params.lat.steps[i]:05.0f} m, step_lon {params.lat.steps[i]:05.0f} m, step_alt {params.alt.steps[i]:03.0f} m, step_off {params.off.steps[i]:04.0f} Hz, step_dft {params.dft.steps[i]:04.0f} Hz,"
                  f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m")

    results.extend(res)
    return IterationResults.found, lat, lon, alt, off, dft, results


def solve(nav_data, satellites, params: CurveFitMethodParameters, init_state: tuple[float, float, float, float, float] | None = None):
    """
    Perform the curve fitting method

    :param nav_data: list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS) | ID
    :param satellites:
    :param params: parameters of solving
    :param init_state: initial state, if None, the function will estimate it
    :return: final lat, lon, alt, off, dft
    """

    # filter curves
    print("SOLVING POSITION")
    detected_curves = find_curves(nav_data,
                                  max_time_gap=params.max_time_gap, min_curve_length=params.min_curve_length)
    print("Detected curves ", len(detected_curves))

    # estimate initial position
    if init_state is None:
        lat_0, lon_0 = estimate_zero_doppler_shift_position(detected_curves, satellites)
        alt_0, off_0, dft_0 = 0, 0, 0
    else:
        lat_0, lon_0, alt_0, off_0, dft_0 = init_state

    # join satellite curves into one array
    curve_array = np.vstack(detected_curves)

    print(f"Initial position: lat {lat_0:05.3f}°, lon {lon_0:05.3f}°, data length {curve_array.shape[0]}")

    # rearrange data
    measured_curve = np.column_stack((curve_array[:, IDX.t],
                                      curve_array[:, IDX.f] - curve_array[:, IDX.fb],
                                      curve_array[:, IDX.fb]))

    # construct arrays
    r = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
             v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
             v_z=curve_array[:, IDX.vz] * unit.km / unit.s)  # for plotting only

    r_sat_arr = np.column_stack((curve_array[:, IDX.x], curve_array[:, IDX.y], curve_array[:, IDX.z]))
    v_sat_arr = np.column_stack((curve_array[:, IDX.vx], curve_array[:, IDX.vy], curve_array[:, IDX.vz]))

    results = list()

    # iterative method
    lat, lon, alt, off, dft = lat_0, lon_0, alt_0, off_0, dft_0
    for _ in range(params.iteration.repeats):
        iter_res, lat, lon, alt, off, dft, results = fit_curve(results, lat, lon, alt, off, dft,
                                                               measured_curve, r_sat_arr, v_sat_arr,
                                                               params.iteration)

    if DEBUG.dump_results:
        dump_data("results", results)
    if DEBUG.plot_results:
        plot_results_of_iterative_position_finding(results, r)

    # grid search method as a followup
    # iter_res, lat, lon, alt, off, dft, results = fit_curve_grid(results, lat, lon, 1500, off, dft,
    #                                                             measured_curve, r_sat_arr, v_sat_arr, params.iteration)
    # if DEBUG.dump_results:
    #     dump_data("results", results)
    # if DEBUG.plot_results:
    #     plot_results_of_iterative_position_finding(results, r)

    if DEBUG.plot_final_curve_fit:
        check_trial_curve(LAT_HOME, LON_HOME, ALT_HOME, 0, 0, measured_curve, r_sat_arr, v_sat_arr, plot=True)

    return lat, lon, alt, off, dft
