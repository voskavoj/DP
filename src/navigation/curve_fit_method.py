from astropy.coordinates import EarthLocation, ITRS
from astropy.time import Time
import astropy.units as unit
import numpy as np
import matplotlib.pyplot as plt

from src.navigation.calculations import m_to_deg_lat, m_to_deg_lon, latlon_distance
from src.satellites.predictions import predict_satellite_positions
from src.navigation.data_processing import nav_data_to_array
from src.navigation.data_processing import NavDataArrayIndices as IDX
from src.config.locations import LOCATIONS
from src.utils.data import dump_data
from src.utils.plots import plot_results_of_iterative_position_finding, plot_analyzed_curve

MIN_CURVE_LEN = 60  # s
MIN_CURVE_DENSITY = 1  # smp/s
MAX_CURVE_VARIANCE = 100  # -
MAX_CURVE_GAP = 5  # s

STEP = 20e3  # km
ITER_LIMIT = 1000
STEP_LIMIT = 10  # m

C = 299792458  # m/s

LON_HOME, LAT_HOME = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1]


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

    # plot_analyzed_curve(curve, dopp_start, dopp_end, curve_duration, curve_density, largest_gap, variance)

    return ((dopp_start > 0 > dopp_end or dopp_start < 0 < dopp_end)
            and curve_duration >= MIN_CURVE_LEN
            and curve_density >= MIN_CURVE_DENSITY
            and largest_gap <= MAX_CURVE_GAP
            and variance <= MAX_CURVE_VARIANCE)


def find_curves(nav_data):
    curves = dict()

    # split by ID and channel
    for time, freq, base_freq, sat_pos, sat_id in nav_data:
        if sat_id not in curves:
            curves[sat_id] = dict()
        if base_freq not in curves[sat_id]:
            curves[sat_id][base_freq] = list()
        curves[sat_id][base_freq].append((time, freq, base_freq, sat_pos, sat_id))

    # to list
    curves = [curves[sat_id][base_freq]
              for sat_id in curves
              for base_freq in curves[sat_id]
              if is_valid_curve(curves[sat_id][base_freq])]

    # todo for actual data
#     # split by gap
#     for curve in curves:
#         prev_time = curve[0][0]
#         for time, freq, base_freq, sat_pos in curve:
#             TODO split by time

    return curves


def solve(nav_data, satellites):
    """

    :param nav_data: list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS) | ID
    :param satellites:
    :return:
    """

    # plot_results_of_iterative_position_finding("results_6", show=True)

    # filter curves
    print("Solving")
    detected_curves = find_curves(nav_data)
    print("Detected curves ", len(detected_curves))
    # for each curve
    for curve in detected_curves:
        sat = satellites[str(curve[0][4])]
        print(sat.name, curve[0][2], len(curve))
        curve_array = nav_data_to_array(curve)

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
        print(Time(pass_time, format="unix").to_value("datetime"), pass_time)

        # find ground position of sat at pass time
        pass_pos_ground = pass_pos.earth_location.geodetic
        print(pass_pos_ground)

        # ITERATIVE ALGORITHM
        iterative_algorithm(curve_array,
                            lat_0=pass_pos_ground.lat.value, lon_0=pass_pos_ground.lon.value, alt_0=0,
                            base_freq=curve[0][IDX.fb])

    plt.show()


def check_trial_curve(lat, lon, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq):
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

    f_d = -1 * rel_vel * base_freq / C
    trial_curve[:, 1] = f_d

    # calculate variance
    sum_of_squares = np.sum((measured_curve[:, 1] - trial_curve[:, 1]) ** 2)

    return sum_of_squares


def fit_curve(results, step, lat_0, lon_0, alt_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq):
    iteration_result = IterationResults.limit

    # iterate
    step_ll = step
    step_alt = 10  # todo
    step_off = 100  # todo
    lat = lat_0
    lon = lon_0
    alt = 0

    # initial SOS
    sos = check_trial_curve(lat, lon, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

    for i in range(ITER_LIMIT):
        # nORTH, sOUTH, wEST, eAST, uP, dOWN, mORE_OFFSET, lESS_OFFSET

        # 1. calculate SOS for each direction
        lat_n = lat + m_to_deg_lat(step_ll * 1, lat)
        sos_n = check_trial_curve(lat_n, lon, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        lat_s = lat + m_to_deg_lat(step_ll * -1, lat)
        sos_s = check_trial_curve(lat_s, lon, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        lon_w = lon + m_to_deg_lon(step_ll * -1, lat)
        sos_w = check_trial_curve(lat, lon_w, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        lon_e = lon + m_to_deg_lon(step_ll * 1, lat)
        sos_e = check_trial_curve(lat, lon_e, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        alt_u = alt + step_alt * 1
        sos_u = check_trial_curve(lat, lon, alt_u, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        alt_d = alt + step_alt * -1
        sos_d = check_trial_curve(lat, lon, alt_d, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        # 2. calculate differences
        sos_n = max(0, sos - sos_n)
        sos_s = max(0, sos - sos_s)
        sos_w = max(0, sos - sos_w)
        sos_e = max(0, sos - sos_e)
        sos_u = max(0, sos - sos_u)
        sos_d = max(0, sos - sos_d)

        diff_lat = sos_n - sos_s
        diff_lon = sos_e - sos_w
        diff_alt = sos_u - sos_d

        # 3. check if current position is the best
        if diff_lat == 0 and diff_lon == 0 and step_ll > STEP_LIMIT:
            step_ll = max(STEP_LIMIT, round(step_ll / 2))

        if diff_lat == 0 and diff_lon == 0 and diff_alt == 0 and step_ll <= STEP_LIMIT:
            results.append([sos, i, lat, lon, alt])
            print(f"Iteration {i + 1:03d}: lat {lat:03.3f}°, lon {lon:02.3f}°, alt {alt:04.0f} m, SOS {sos:09.0f}, "
                  f"step_LL {step_ll:05.0f} m, "
                  f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m, "
                  f"ITERATION END")
            if i == 0:
                iteration_result = IterationResults.first
            else:
                iteration_result = IterationResults.found
            break

        # 4. calculate ratios
        tot_diff = abs(diff_lat) + abs(diff_lon) + abs(diff_alt)
        if tot_diff > 0:
            ratio_lat = diff_lat / tot_diff
            ratio_lon = diff_lon / tot_diff
            ratio_alt = diff_alt / tot_diff
        else:
            ratio_lat, ratio_lon, ratio_alt = 0, 0, 0

        # 5. move trial position and calculate new SOS
        lat = lat + m_to_deg_lat(step_ll * ratio_lat, lat)
        lon = lon + m_to_deg_lon(step_ll * ratio_lon, lat)
        alt = alt + step_alt * ratio_alt

        sos = check_trial_curve(lat, lon, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        # 6. calculate new step
        # step_ll = max(step_ll, 3 * sos / 1000)  # adaptive step
        # step_ll = max(min(STEP, step_ll), STEP_LIMIT)  # constrain step

        # 7. log results
        results.append([sos, i, lat, lon, alt])
        print(f"Iteration {i + 1:03d}: lat {lat:05.3f}°, lon {lon:05.3f}°, alt {alt:04.0f} m, SOS {sos:09.0f}, "
              f"step_LL {step_ll:05.0f} m, "
              f"dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):07.0f} m, "
              f"rlat {ratio_lat: 5.2f}, rlon {ratio_lon: 5.2f}, ralt {ratio_alt: 5.2f}")

    return iteration_result, lat, lon, alt, results


def iterative_algorithm(curve_array, lat_0, lon_0, alt_0, base_freq):

    measured_curve = np.column_stack((curve_array[:, IDX.t], curve_array[:, IDX.f] - curve_array[:, IDX.fb]))

    time_arr = Time(curve_array[:, IDX.t], format="unix")
    r = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
             v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
             v_z=curve_array[:, IDX.vz] * unit.km / unit.s,
             obstime=time_arr)
    r_sat_arr = r.cartesian.without_differentials()
    v_sat_arr = r.velocity

    results = list()
    step = STEP
    lat, lon, alt = lat_0, lon_0, alt_0

    while step > STEP_LIMIT:
        print(f"Step {step}")
        iter_res, lat, lon, alt, results = fit_curve(results, step, lat, lon, alt,
                                                     measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
        # step = round(step / 2)

        if iter_res == IterationResults.found:
            break

        # if we hit iteration limit, it makes no sense to continue
        if iter_res == IterationResults.limit:
            print("Iteration limit reached")
            break

    dump_data("results", results)
    plot_results_of_iterative_position_finding(results, r)
    plt.show()

    final_lat, final_lon = min(results, key=lambda x: x[0])[2], min(results, key=lambda x: x[0])[3]
    print(f"Position error: {latlon_distance(LAT_HOME, final_lat, LON_HOME, final_lon):.1f} m")
