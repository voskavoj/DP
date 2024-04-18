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
from src.config.locations import LOCATIONS
from src.utils.data import dump_data
from src.utils.plots import plot_results_of_iterative_position_finding, plot_analyzed_curve, \
    plot_measured_vs_trial_curve

# Curve validation parameters (0 to disable)
MIN_CURVE_LEN = 10  # samples
MAX_CURVE_GAP = 60  # s

# Constraints
ALT_LIMIT = 3000  # m
STEP_LIMIT = 10  # m

# Iterative algorithm parameters
ITER_INIT_STEP_LL = 100e3  # km
ITER_INIT_STEP_ALT = 0  # m
ITER_INIT_STEP_OFF = 3000  # Hz
ITER_LIMIT = 500
ITER_INIT_STEP_DFT = 100  # mHz/s

# Grid search parameters
GRID_SEARCH_DEPTH = 5
GRID_SIZE = 10
GRID_INIT_STEP_LL = 100e3  # km
GRID_INIT_STEP_ALT = 100  # m
GRID_INIT_STEP_OFF = 3000  # Hz

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

        :return: lat, lon; or 0, 0 if no valid curves were found
    """

    estimated_init_locations = list()
    for curve_array in detected_curves:
        sat = satellites[str(int(curve_array[0, IDX.sat_id]))]
        print(sat.name, curve_array[0, IDX.fb], curve_array.shape[0])

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
    f_drift = (trial_curve[:, 0] - np.min(trial_curve[:, 0])) * (dft / 1000)  # todo make it so it's in base units - individual limits!

    # Construct trial curve
    trial_curve[:, 1] = f_d + off + f_drift

    # Calculate variance
    sum_of_squares = np.sum((measured_curve[:, 1] - trial_curve[:, 1]) ** 2)

    if plot:
        plot_measured_vs_trial_curve(measured_curve, trial_curve, lat, lon, alt, off)

    return sum_of_squares


def _adjust_step(step, diff, force=False):
    if (force or diff == 0) and step > STEP_LIMIT:
        step = max(STEP_LIMIT, round(step * 0.75))
    return step


class RatioCycleDetector:
    def __init__(self):
        self.last_ratios = None
        self.second_to_last_ratios = None

    def update(self, ratios):
        self.second_to_last_ratios = self.last_ratios
        self.last_ratios = ratios

    def is_cycle(self):
        if self.last_ratios is None or self.second_to_last_ratios is None:
            return False

        cycle = [self.last_ratios[i] + self.second_to_last_ratios[i] for i in range(len(self.last_ratios))] == [0, 0, 0, 0]

        if cycle:
            self.last_ratios, self.second_to_last_ratios = None, None

        return cycle


def fit_curve(results, lat_0, lon_0, alt_0, off_0, dft_0, measured_curve, r_sat_arr, v_sat_arr):
    iteration_result = IterationResults.limit
    cycle_detector = RatioCycleDetector()

    # iterate
    step_lat, step_lon = ITER_INIT_STEP_LL, ITER_INIT_STEP_LL
    step_alt = ITER_INIT_STEP_ALT
    step_off = ITER_INIT_STEP_OFF
    step_dft = ITER_INIT_STEP_DFT
    lat = lat_0
    lon = lon_0
    alt = alt_0
    off = off_0
    dft = dft_0

    # initial SOS
    sos = check_trial_curve(lat, lon, alt, off, 0, measured_curve, r_sat_arr, v_sat_arr)
    results.append([sos, -1, lat, lon, alt])

    for i in range(ITER_LIMIT):
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

        if (alt == 0 and diff_alt < 0) or (alt == ALT_LIMIT and diff_alt > 0):
            diff_alt = 0

        step_lat = _adjust_step(step_lat, diff_lat)
        step_lon = _adjust_step(step_lon, diff_lon)
        step_alt = _adjust_step(step_alt, diff_alt)
        step_off = _adjust_step(step_off, diff_off)
        step_dft = _adjust_step(step_dft, diff_dft)

        if (diff_lat == 0 and diff_lon == 0 and diff_off == 0 and diff_alt == 0 and diff_dft == 0
                and step_lat <= STEP_LIMIT and step_lon <= STEP_LIMIT):

            results.append([sos, i, lat, lon, alt])
            print(f"Iteration {i + 1:03d}: lat {lat:03.3f}°, lon {lon:02.3f}°, alt {alt:04.0f} m, SOS {sos:09.0f}, "
                  f"step_lat {step_lat:05.0f} m, step_lon {step_lon:05.0f} m, "
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
        lat = lat + m_to_deg_lat(step_lat * ratio_lat, lat)
        lon = lon + m_to_deg_lon(step_lon * ratio_lon, lat)
        alt = min(ALT_LIMIT, max(0, alt + step_alt * ratio_alt))
        off = off + step_off * ratio_off
        dft = dft + step_dft * ratio_dft

        # check if the algorithm is not stuck in a cycle
        cycle_detector.update([ratio_lat, ratio_lon, ratio_alt, ratio_off])
        if cycle_detector.is_cycle():  # cycle detection
            step_lat = _adjust_step(step_lat, diff_lat, force=ratio_lat != 0)
            step_lon = _adjust_step(step_lon, diff_lon, force=ratio_lon != 0)
            step_alt = _adjust_step(step_alt, diff_alt, force=ratio_alt != 0)
            step_off = _adjust_step(step_off, diff_off, force=ratio_off != 0)
            # step_dft = _adjust_step(step_dft, diff_dft, force=ratio_dft != 0)

            lat = lat + m_to_deg_lat(step_lat * ratio_lat, lat)
            lon = lon + m_to_deg_lon(step_lon * ratio_lon, lat)
            alt = min(ALT_LIMIT, max(0, alt + step_alt * ratio_alt))
            off = off + step_off * ratio_off
            dft = dft + step_dft * ratio_dft

        sos = check_trial_curve(lat, lon, alt, off, dft, measured_curve, r_sat_arr, v_sat_arr)

        # 6. log results
        results.append([sos, i, lat, lon, alt])
        print(f"Iteration {i + 1:03d}: lat {lat:05.3f}°, lon {lon:05.3f}°, alt {alt:04.0f} m, off {off:05.0f}, dft {dft:04.0f}, SOS {sos:09.0f}, "
              f"step_lat {step_lat:05.0f} m, step_lon {step_lon:05.0f} m, step_alt {step_alt:03.0f} m, step_off {step_off:04.0f} Hz, step_dft {step_dft:04.0f} mHz/s,"
              f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m, "
              f"rlat {ratio_lat: 2.0f}, rlon {ratio_lon: 2.0f}, ralt {ratio_alt: 2.0f}, roff {ratio_off: 2.0f}, rdft {ratio_dft: 2.0f}{', Cycle' if cycle_detector.is_cycle() else ''}")

    # check_trial_curve(lat, lon, alt, off, 0, measured_curve, r_sat_arr, v_sat_arr, plot=True)
    return iteration_result, lat, lon, alt, off, dft, results


def fit_curve_grid(results, lat_0, lon_0, alt_0, off_0, dft_0, measured_curve, r_sat_arr, v_sat_arr):
    iteration_result = IterationResults.limit

    # iterate
    step_lat, step_lon = [GRID_INIT_STEP_LL] * 2
    step_alt = GRID_INIT_STEP_ALT
    step_off = GRID_INIT_STEP_OFF
    lat = lat_0
    lon = lon_0
    alt = alt_0
    off = off_0
    dft = dft_0

    res = list()

    for i in range(GRID_SEARCH_DEPTH):
        # grid: 0 = lat, 1 = lon, 2 = alt, 3 = off
        lats = [lat - m_to_deg_lat(step_lat * k, lat) for k in range(-GRID_SIZE, GRID_SIZE + 1)] if step_lat >= STEP_LIMIT else [lat]
        lons = [lon - m_to_deg_lat(step_lon * k, lat) for k in range(-GRID_SIZE, GRID_SIZE + 1)] if step_lon >= STEP_LIMIT else [lon]
        # alts = [alt - step_alt * k for k in range(-GRID_SIZE, GRID_SIZE + 1)] if step_alt >= STEP_LIMIT else [alt]  # todo way too slow
        alts = [0]
        offs = [off - step_off * k for k in range(-GRID_SIZE, GRID_SIZE + 1)] if step_off >= STEP_LIMIT else [off]

        for lat, lon, alt, off in itertools.product(lats, lons, alts, offs):
            sos = check_trial_curve(lat, lon, alt, off, 0, measured_curve, r_sat_arr, v_sat_arr)
            # print(sos, lat, lon, alt, off)
            res.append([sos, 0, lat, lon, alt, off])

        sos, _, lat, lon, alt, off = min(res, key=lambda x: x[0])
        res.append([sos, 0, lat, lon, alt, off])
        print(f"Grid {i + 1:03d}: lat {lat:05.3f}°, lon {lon:05.3f}°, alt {alt:04.0f} m, off {off:05.0f} SOS {sos:09.0f}, "
              f"step_lat {step_lat:05.0f} m, step_lon {step_lon:05.0f} m, step_alt {step_alt:03.0f} m, step_off {step_off:04.0f} Hz,"
              f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m")

        step_lat, step_lon, step_alt, step_off = round(step_lat / GRID_SIZE), round(step_lon / GRID_SIZE), round(step_alt / GRID_SIZE), round(step_off / GRID_SIZE)

    return IterationResults.found, lat, lon, alt, off, dft, res


def iterative_algorithm(curve_array, lat_0, lon_0, alt_0, off_0, dft_0):

    measured_curve = np.column_stack((curve_array[:, IDX.t], curve_array[:, IDX.f] - curve_array[:, IDX.fb], curve_array[:, IDX.fb]))

    r = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
             v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
             v_z=curve_array[:, IDX.vz] * unit.km / unit.s)

    r_sat_arr = np.column_stack((curve_array[:, IDX.x], curve_array[:, IDX.y], curve_array[:, IDX.z]))
    v_sat_arr = np.column_stack((curve_array[:, IDX.vx], curve_array[:, IDX.vy], curve_array[:, IDX.vz]))

    results = list()

    iter_res, lat, lon, alt, off, dft, results = fit_curve(results, lat_0, lon_0, alt_0, off_0, dft_0,
                                                           measured_curve, r_sat_arr, v_sat_arr)
    dump_data("results", results)
    plot_results_of_iterative_position_finding(results, r)

    # iter_res, lat, lon, alt, off, dft, results = fit_curve_grid(results, lat_0, lon_0, alt_0, off_0, dft_0,
    #                                                             measured_curve, r_sat_arr, v_sat_arr)
    # dump_data("results", results)
    # plot_results_of_iterative_position_finding(results, r)


def solve(nav_data, satellites):
    """

    :param nav_data: list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS) | ID
    :param satellites:
    :return:
    """

    # filter curves
    print("Solving")
    detected_curves = find_curves(nav_data, max_time_gap=MAX_CURVE_GAP, min_curve_length=MIN_CURVE_LEN)
    print("Detected curves ", len(detected_curves))

    # estimate initial position
    lat_0, lon_0 = estimate_zero_doppler_shift_position(detected_curves, satellites)
    curve_array = np.vstack(detected_curves)
    print(curve_array.shape)

    print(f"Initial position: lat {lat_0:05.3f}°, lon {lon_0:05.3f}°")
    iterative_algorithm(curve_array,
                        lat_0=lat_0, lon_0=lon_0, alt_0=0, off_0=0, dft_0=0)
