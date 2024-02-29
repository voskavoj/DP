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
ITER_LIMIT = 150
STEP_LIMIT = 10  # m

C = 299792458  # m/s

LON_HOME, LAT_HOME = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1]


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

    plot_analyzed_curve(curve, dopp_start, dopp_end, curve_duration, curve_density, largest_gap, variance)

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
        try:
            iterative_algorithm(curve_array,
                                lat_0=pass_pos_ground.lat.value, lon_0=pass_pos_ground.lon.value, alt_0=0,
                                base_freq=curve[0][IDX.fb])
        except KeyboardInterrupt:
            print("Interrupted")

    plt.show()


def check_trial_curve(lat, lon, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq):
    curve_len = measured_curve.shape[0]

    trial_curve = np.empty((curve_len, 2))
    trial_curve[:, 0] = measured_curve[:, 0]  # times are the same

    r_user_arr = (EarthLocation.from_geodetic([lon] * curve_len, [lat] * curve_len, [alt] * curve_len)
                  .get_itrs(obstime=time_arr).cartesian.without_differentials())

    # calculate doppler curve
    for k in range(curve_len):
        vs, rs, ru = v_sat_arr[k], r_sat_arr[k], r_user_arr[k]
        vs = np.array([vs.d_x.value, vs.d_y.value, vs.d_z.value]) * 1000
        rs = np.array([rs.x.value, rs.y.value, rs.z.value]) * 1000
        ru = np.array([ru.x.value, ru.y.value, ru.z.value]) * 1
        rel_vel = np.dot(vs, (rs - ru) / np.linalg.norm(rs - ru))

        f_d = -1 * rel_vel * base_freq / C
        trial_curve[k, 1] = f_d

    # calculate variance
    sum_of_squares = np.sum((measured_curve[:, 1] - trial_curve[:, 1]) ** 2)

    return sum_of_squares


def move_latlon(step, lat, lon, north=False, south=False, west=False, east=False):
    if (north and south) or (west and east) or not any([north, south, west, east]):
        raise ValueError(f"Invalid direction combination (NSWE): {north, south, west, east}")

    if north:
        lat += m_to_deg_lat(step, lat)
    if south:
        lat -= m_to_deg_lat(step, lat)
    if west:
        lon -= m_to_deg_lon(step, lat)
    if east:
        lon += m_to_deg_lon(step, lat)

    return lat, lon


def fit_curve(results, step, lat_0, lon_0, alt_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq):

    # init
    go_north, go_south, go_west, go_east = False, False, False, False
    # first, get SOS for current location
    sos_0 = check_trial_curve(lat_0, lon_0, alt_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
    # second, go north, south, west and east to see the direction
    lat, lon = move_latlon(step, lat_0, lon_0, north=True)
    sos_n = check_trial_curve(lat, lon, alt_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
    if sos_n < sos_0:
        go_north = True
    else:
        lat, lon = move_latlon(step, lat_0, lon_0, south=True)
        sos_s = check_trial_curve(lat, lon, alt_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
        if sos_s < sos_0:
            go_south = True

    lat, lon = move_latlon(step, lat_0, lon_0, west=True)
    sos_w = check_trial_curve(lat, lon, alt_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
    if sos_w < sos_0:
        go_west = True
    else:
        lat, lon = move_latlon(step, lat_0, lon_0, east=True)
        sos_e = check_trial_curve(lat, lon, alt_0, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
        if sos_e < sos_0:
            go_east = True

    print(f"Go N {go_north}, S {go_south}, W {go_west}, E {go_east}")

    if not any([go_north, go_south, go_west, go_east]):
        print("Original location is the best")
        return lat_0, lon_0, alt_0, results

    # iterate
    lat = lat_0
    lon = lon_0
    alt = 0
    sos = sos_0

    for i in range(ITER_LIMIT):
        # check north-south
        if go_north or go_south:
            lat_ns, lon_ns = move_latlon(step, lat, lon, north=go_north, south=go_south)
            sos_ns = check_trial_curve(lat_ns, lon_ns, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
            if sos_ns < sos:  # we should move in ns direction
                lat = lat_ns
            else:
                go_north, go_south = False, False

        # check west-east
        if go_west or go_east:
            lat_we, lon_we = move_latlon(step, lat, lon, west=go_west, east=go_east)
            sos_we = check_trial_curve(lat_we, lon_we, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
            if sos_we < sos:
                lon = lon_we
            else:
                go_west, go_east = False, False

        if not any([go_north, go_south, go_west, go_east]):
            print("Found position!")
            break

        sos = check_trial_curve(lat, lon, alt, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)

        results.append([sos, i, lat, lon, alt])
        print(f"Iteration {i + 1:03d}: lat {lat:03.3f}, lon {lon:02.3f}, alt {alt:04.0f}, SOS {sos:09.0f}, "
              f"{latlon_distance(LAT_HOME, lat, LON_HOME, lon):.0f} m"
              f"{'N' if go_north else ' '}{'S' if go_south else ' '}{'W' if go_west else ' '}{'E' if go_east else ' '}")

    return lat, lon, alt, results


def iterative_algorithm(curve_array, lat_0, lon_0, alt_0, base_freq):

    measured_curve = np.column_stack((curve_array[:, IDX.t], curve_array[:, IDX.f] - curve_array[:, IDX.fb]))

    time_arr = Time(curve_array[:, IDX.t], format="unix")
    r = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
             v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
             v_z=curve_array[:, IDX.vz] * unit.km / unit.s,
             obstime=time_arr)
    r_sat_arr = r.cartesian.without_differentials()
    v_sat_arr = r.velocity

    plt.figure()
    plt.plot(measured_curve[:, 0], measured_curve[:, 1], "-k", linewidth=0.5)

    results = list()
    step = STEP
    lat, lon, alt = lat_0, lon_0, alt_0

    while step > STEP_LIMIT:
        print(f"Step {step}")
        lat, lon, alt, results = fit_curve(results, step, lat, lon, alt,
                                           measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq)
        step = int(step / 2)

    dump_data("results", results)
    plot_results_of_iterative_position_finding(results, r)

    final_lat, final_lon = min(results, key=lambda x: x[0])[2], min(results, key=lambda x: x[0])[3]
    print(f"Position error: {latlon_distance(LAT_HOME, final_lat, LON_HOME, final_lon):.1f} m")
