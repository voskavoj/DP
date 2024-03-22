import itertools

from astropy.coordinates import EarthLocation, ITRS
from astropy.time import Time
import astropy.units as unit
import numpy as np
import matplotlib.pyplot as plt

from src.navigation.calculations import m_to_deg_lat, m_to_deg_lon, latlon_distance
from src.satellites.predictions import predict_satellite_positions
from src.navigation.data_processing import nav_data_to_array, find_curves
from src.navigation.data_processing import NavDataArrayIndices as IDX
from src.config.locations import LOCATIONS
from src.utils.data import dump_data
from src.utils.plots import plot_results_of_iterative_position_finding, plot_analyzed_curve

# Curve validation parameters (False to disable)
TREND = True
MIN_CURVE_LEN = 30  # s
MIN_CURVE_DENSITY = 1  # smp/s
MAX_CURVE_VARIANCE = 250  # -
MAX_CURVE_GAP = False  # 60  # s

# Constraints
ALT_LIMIT = 3000  # m
STEP_LIMIT = 10  # m

# Iterative algorithm parameters
ITER_INIT_STEP_LL = 100e3  # km
ITER_INIT_STEP_ALT = 100  # m
ITER_INIT_STEP_OFF = 3000  # Hz
ITER_LIMIT = 500

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


def is_valid_curve(curve):
    """
    Determine if a curve is valid
    Criteria:   min duration (MIN_CURVE_LEN)
                max variance (MAX_CURVE_VARIANCE)
                min density (MIN_CURVE_DENSITY)
                max gap (MAX_CURVE_GAP)
                trend (doppler shift must change sign)
                must be longer than MIN_CURVE_LEN
    :param curve:
    :return: bool
    """

    try:
        dopp_start = curve[0][1] - curve[0][2]
        dopp_end = curve[-1][1] - curve[-1][2]
        curve_duration = (curve[-1][0] - curve[0][0]).sec
        curve_density = 0 if not curve_duration else len(curve) / curve_duration
        largest_gap = max([(curve[i][0] - curve[i - 1][0]).sec for i in range(1, len(curve))])
        variance = sum([abs(curve[i - 1][1] - curve[i][1]) for i in range(1, len(curve))]) / len(curve)
    except (ValueError, RuntimeError):
        return False

    valid = bool(
        (TREND is False or (dopp_start > 0 > dopp_end or dopp_start < 0 < dopp_end))
        and (MIN_CURVE_LEN is False or curve_duration >= MIN_CURVE_LEN)
        and (MIN_CURVE_DENSITY is False or curve_density >= MIN_CURVE_DENSITY)
        and (MAX_CURVE_GAP is False or largest_gap <= MAX_CURVE_GAP)
        and (MAX_CURVE_VARIANCE is False or variance <= MAX_CURVE_VARIANCE)
    )

    plot_analyzed_curve(curve, dopp_start, dopp_end, curve_duration, curve_density, largest_gap, variance, ok=valid)

    return valid


def solve(nav_data, satellites):
    """

    :param nav_data: list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS) | ID
    :param satellites:
    :return:
    """

    # filter curves
    print("Solving")
    detected_curves = find_curves(nav_data)
    print("Detected curves ", len(detected_curves))

    init_locations = list()
    measured_curves = list()

    for curve_array in detected_curves:
        sat = satellites[str(int(curve_array[0, IDX.sat_id]))]
        print(sat.name, curve_array[0, IDX.fb], curve_array.shape[0])

        # find zero doppler shift
        dopp_shift_arr = curve_array[:, IDX.f] - curve_array[:, IDX.fb]
        if dopp_shift_arr[0] > dopp_shift_arr[-1]:  # for interpolation the doppler shift must be increasing
            dopp_shift_arr = dopp_shift_arr[::-1]
        pass_time = np.interp(0, xp=dopp_shift_arr, fp=curve_array[:, IDX.t], left=0, right=0)
        if pass_time == 0:
            print("SKIPPED")
            # continue

        times = pass_time + np.array([-30, 0, 30])
        pred_pos = predict_satellite_positions([sat.satrec], Time(times, format="unix"))[0, :]
        pass_pos = pred_pos[1]

        # find ground position of sat at pass time
        pass_pos_ground = pass_pos.earth_location.geodetic

        init_locations.append((pass_pos_ground.lat.value, pass_pos_ground.lon.value, 0))
        measured_curves.append(curve_array)

    lat_0, lon_0, alt_0 = 0, 0, 0
    for lat, lon, alt in init_locations:
        lat_0 += lat
        lon_0 += lon
        alt_0 += alt
    lat_0 /= len(init_locations)
    lon_0 /= len(init_locations)
    alt_0 /= len(init_locations)

    curve_array = np.vstack([curve for curve in measured_curves if curve[0, IDX.fb] == 1626270800])

    print(curve_array.shape)

    print(f"Initial position: lat {lat_0:05.3f}°, lon {lon_0:05.3f}°, alt {alt_0:04.0f} m")
    iterative_algorithm(curve_array,
                        lat_0=lat_0, lon_0=lon_0, alt_0=alt_0, off_0=0,
                        base_freq=1626270800.0)  # todo different base freqs


def check_trial_curve(lat, lon, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq):
    curve_len = measured_curve.shape[0]

    trial_curve = np.empty((curve_len, 2))
    trial_curve[:, 0] = measured_curve[:, 0]  # times are the same
    r_user_arr = (EarthLocation.from_geodetic(lon, lat, alt)
                  .get_itrs(obstime=time_arr).cartesian.without_differentials())

    vs, rs, ru = v_sat_arr, r_sat_arr, r_user_arr
    vs = np.array([vs.d_x.value, vs.d_y.value, vs.d_z.value]) * 1000
    rs = np.array([rs.x.value, rs.y.value, rs.z.value]) * 1000
    ru = np.array([ru.x.value, ru.y.value, ru.z.value]) * 1
    # Ro_dot = (V_s - V_u) * (r_s - r_u) / ||r_s - r_u||
    rel_vel = np.sum(vs * (rs - ru) / np.linalg.norm(rs - ru, axis=0), axis=0)

    f_b = measured_curve[:, 2]
    f_d = -1 * rel_vel * f_b / C
    trial_curve[:, 1] = f_d + off

    # calculate variance
    sum_of_squares = np.sum((measured_curve[:, 1] - trial_curve[:, 1]) ** 2)

    return sum_of_squares


def fit_curve(results, lat_0, lon_0, alt_0, off_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq):
    iteration_result = IterationResults.limit

    # iterate
    step_lat, step_lon = ITER_INIT_STEP_LL, ITER_INIT_STEP_LL
    step_alt = ITER_INIT_STEP_ALT
    step_off = ITER_INIT_STEP_OFF
    lat = lat_0
    lon = lon_0
    alt = alt_0
    off = off_0

    # initial SOS
    sos = check_trial_curve(lat, lon, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

    for i in range(ITER_LIMIT):
        # nORTH, sOUTH, wEST, eAST, uP, dOWN, mORE_OFFSET, lESS_OFFSET

        # 1. calculate SOS for each direction
        lat_n = lat + m_to_deg_lat(step_lat * 1, lat)
        sos_n = check_trial_curve(lat_n, lon, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        lat_s = lat + m_to_deg_lat(step_lat * -1, lat)
        sos_s = check_trial_curve(lat_s, lon, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        lon_w = lon + m_to_deg_lon(step_lon * -1, lat)
        sos_w = check_trial_curve(lat, lon_w, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        lon_e = lon + m_to_deg_lon(step_lon * 1, lat)
        sos_e = check_trial_curve(lat, lon_e, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        alt_u = alt + step_alt * 1
        sos_u = check_trial_curve(lat, lon, alt_u, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        alt_d = alt + step_alt * -1
        sos_d = check_trial_curve(lat, lon, alt_d, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
        
        off_m = off + step_off * 1
        sos_m = check_trial_curve(lat, lon, alt, off_m, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        off_l = off + step_off * -1
        sos_l = check_trial_curve(lat, lon, alt, off_l, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        # 2. calculate differences
        sos_n = max(0, sos - sos_n)
        sos_s = max(0, sos - sos_s)
        sos_w = max(0, sos - sos_w)
        sos_e = max(0, sos - sos_e)
        sos_u = max(0, sos - sos_u)
        sos_d = max(0, sos - sos_d)
        sos_m = max(0, sos - sos_m)
        sos_l = max(0, sos - sos_l)

        diff_lat = sos_n - sos_s
        diff_lon = sos_e - sos_w
        diff_alt = sos_u - sos_d
        diff_off = sos_m - sos_l

        # 3. check if current position is the best
        if (alt == 0 and diff_alt < 0) or (alt == ALT_LIMIT and diff_alt > 0):
            diff_alt = 0

        if diff_lat == 0 and step_lat > STEP_LIMIT:
            step_lat = max(STEP_LIMIT, round(step_lat * 0.75))
        if diff_lon == 0 and step_lon > STEP_LIMIT:
            step_lon = max(STEP_LIMIT, round(step_lon * 0.75))
        if diff_alt == 0 and step_alt > STEP_LIMIT:
            step_alt = max(STEP_LIMIT, round(step_alt * 0.75))
        if diff_off == 0 and step_off > STEP_LIMIT:
            step_off = max(STEP_LIMIT, round(step_off * 0.75))

        if (diff_lat == 0 and diff_lon == 0 and diff_off == 0 and diff_alt == 0
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

        # 5. move trial position and calculate new SOS
        lat = lat + m_to_deg_lat(step_lat * ratio_lat, lat)
        lon = lon + m_to_deg_lon(step_lon * ratio_lon, lat)
        alt = min(ALT_LIMIT, max(0, alt + step_alt * ratio_alt))
        off = off + step_off * ratio_off

        sos = check_trial_curve(lat, lon, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        # 6. log results
        results.append([sos, i, lat, lon, alt])
        print(f"Iteration {i + 1:03d}: lat {lat:05.3f}°, lon {lon:05.3f}°, alt {alt:04.0f} m, off {off:05.0f} SOS {sos:09.0f}, "
              f"step_lat {step_lat:05.0f} m, step_lon {step_lon:05.0f} m, step_alt {step_alt:03.0f} m, step_off {step_off:05.0f} Hz,"
              f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m, "
              f"rlat {ratio_lat: 2.0f}, rlon {ratio_lon: 2.0f}, ralt {ratio_alt: 2.0f}, roff {ratio_off: 2.0f}")

    return iteration_result, lat, lon, alt, results


def fit_curve_grid(results, lat_0, lon_0, alt_0, off_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq):
    iteration_result = IterationResults.limit

    # iterate
    step_lat, step_lon = [GRID_INIT_STEP_LL] * 2
    step_alt = GRID_INIT_STEP_ALT
    step_off = GRID_INIT_STEP_OFF
    lat = lat_0
    lon = lon_0
    alt = alt_0
    off = off_0

    res = list()

    for i in range(GRID_SEARCH_DEPTH):
        # grid: 0 = lat, 1 = lon, 2 = alt, 3 = off
        lats = [lat - m_to_deg_lat(step_lat * k, lat) for k in range(-GRID_SIZE, GRID_SIZE + 1)]
        lons = [lon - m_to_deg_lat(step_lon * k, lat) for k in range(-GRID_SIZE, GRID_SIZE + 1)]
        alts = [alt - step_alt * k for k in range(-GRID_SIZE, GRID_SIZE + 1)]
        offs = [off - step_off * k for k in range(-GRID_SIZE, GRID_SIZE + 1)]

        for lat, lon, alt, off in itertools.product(lats, lons, [0], offs):
            sos = check_trial_curve(lat, lon, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
            # print(sos, lat, lon, alt, off)
            res.append([sos, 0, lat, lon, alt, off])

        sos, _, lat, lon, alt, off = min(res, key=lambda x: x[0])
        res.append([sos, 0, lat, lon, alt, off])
        print(f"Grid {i + 1:03d}: lat {lat:05.3f}°, lon {lon:05.3f}°, alt {alt:04.0f} m, off {off:05.0f} SOS {sos:09.0f}, "
              f"step_lat {step_lat:05.0f} m, step_lon {step_lon:05.0f} m, step_alt {step_alt:03.0f} m, step_off {step_off:05.0f} Hz,"
              f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m")

        step_lat, step_lon, step_alt, step_off = round(step_lat / GRID_SIZE), round(step_lon / GRID_SIZE), round(step_alt / GRID_SIZE), round(step_off / GRID_SIZE)

    return IterationResults.found, lat, lon, alt, res


def iterative_algorithm(curve_array, lat_0, lon_0, alt_0, off_0, base_freq):

    measured_curve = np.column_stack((curve_array[:, IDX.t], curve_array[:, IDX.f] - curve_array[:, IDX.fb], curve_array[:, IDX.fb]))

    time_arr = Time(curve_array[:, IDX.t], format="unix")
    r = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
             v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
             v_z=curve_array[:, IDX.vz] * unit.km / unit.s,
             obstime=time_arr)
    r_sat_arr = r.cartesian.without_differentials()
    v_sat_arr = r.velocity

    results = list()

    iter_res, lat, lon, alt, results = fit_curve(results, lat_0, lon_0, alt_0, off_0,
                                                 measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
    dump_data("results", results)
    plot_results_of_iterative_position_finding(results, r)

    iter_res, lat, lon, alt, results = fit_curve_grid(results, lat_0, lon_0, alt_0, off_0,
                                                      measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
    dump_data("results", results)
    plot_results_of_iterative_position_finding(results, r)

    print(f"Position error: {latlon_distance(LAT_HOME, lat, LON_HOME, lon, ALT_HOME, alt):.1f} m")
