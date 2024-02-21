import random

import numpy as np
from astropy.coordinates import EarthLocation, ITRS, CartesianRepresentation, CartesianDifferential
from astropy.time import Time, TimeDelta
import astropy.units as unit
import matplotlib.pyplot as plt
import math

from src.radio.iridium_offline_radio import IridiumOfflineRadio
from src.satellites.download_tle import download_tles
from src.navigation.data_processing import process_received_frames

from src.radio.iridium_channels import map_sat_id_to_tle_id
from src.satellites.predictions import predict_satellite_positions
from src.navigation.data_processing import nav_data_to_array
from src.navigation.data_processing import NavDataArrayIndices as IDX


FIG_IDX = 100

MIN_CURVE_LEN = 60  # s
MIN_CURVE_DENSITY = 1  # smp/s
MAX_CURVE_VARIANCE = 100  # -
MAX_CURVE_GAP = 5  # s

STEP = 10  # km
ITER_LIMIT = 100

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

    # print
    curve_array = nav_data_to_array(curve)
    plt.figure()
    plt.title(f"T:{(dopp_start > 0 > dopp_end or dopp_start < 0 < dopp_end)}, "
              f"L:{curve_duration:.1f}, "
              f"D:{curve_density:.3f} \n"
              f"G{largest_gap:.2f}, "
              f"V{variance:.3f}\n"
              f"{((dopp_start > 0 > dopp_end or dopp_start < 0 < dopp_end)
                 and curve_duration >= MIN_CURVE_LEN
                 and curve_density >= MIN_CURVE_DENSITY
                 and largest_gap <= MAX_CURVE_GAP
                 and variance <= MAX_CURVE_VARIANCE)}")
    plt.plot(curve_array[:, 0], curve_array[:, 1] - curve_array[:, 2], ".")
    global FIG_IDX
    plt.savefig(f"f{FIG_IDX}.png")
    plt.close()
    FIG_IDX += 1

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
    detected_curves = find_curves(nav_data)
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

        # plt.figure()
        # plt.plot(curve_array[:, 0], dopp_shift_arr, ".")
        # plt.plot(pass_time, 0, "r.")

        # find ground position of sat at pass time
        pass_pos_ground = pass_pos.earth_location.geodetic
        print(pass_pos_ground)

        # calculate satellite asimuth (of travel)
        pred_pos_ground = pred_pos.earth_location.geodetic
        d_lon = pred_pos_ground.lon[0] - pred_pos_ground.lon[-1]
        d_lat = pred_pos_ground.lat[0] - pred_pos_ground.lat[-1]
        sat_az = np.arctan2(d_lon, d_lat).value  # todo check
        print(sat_az, np.rad2deg(sat_az))

        sat_az += np.pi / 2  # perpendicular to travel
        print(sat_az, np.rad2deg(sat_az))
        # generate coarse list of possible locations

        # ITERATIVE ALGORITHM
        iterative_algorithm(curve, curve_array, lat_0=pass_pos_ground.lat.value, lon_0=pass_pos_ground.lon.value, az=sat_az)

        # for each location generate a doppler curve based on the timestamps of the original curves
        # find the best match
    # plt.show()


def geodetic_distance(angle_dist, coord):
    """
    Approximate conversion from distance in angle to distance in km based on latitude/longitude

    https://en.wikipedia.org/wiki/Latitude#Meridian_distance_on_the_ellipsoid

    :param angle_dist:
    :param coord:
    :return:
    """
    EARTH_RADIUS = 6371  # km
    m_per_deg_lat = 111132.954 - 559.822 * np.cos(2 * latMid) + 1.175 * np.cos(4 * latMid)
    m_per_deg_lon = 111132.954 * np.cos(latMid)


def m_to_deg_lat(m, lat, deg=True):
    """
    Approximate conversion from distance in km to distance in latitude
    https://en.wikipedia.org/wiki/Latitude#Meridian_distance_on_the_ellipsoid

    :param m: meters
    :param lat: current latitude (deg)
    :param deg: input is in degrees (True) or radians (False)
    :return: degrees latitude
    """
    if deg:
        lat = np.deg2rad(lat)

    m_per_deg_lat = 111132.954 - 559.822 * np.cos(2 * lat) + 1.175 * np.cos(4 * lat)
    return m / m_per_deg_lat


def m_to_deg_lon(m, lat, deg=True):
    """
    Approximate conversion from distance in km to distance in longitude
    https://en.wikipedia.org/wiki/Latitude#Meridian_distance_on_the_ellipsoid

    :param m: meters
    :param lat: current latitude (deg)
    :param deg: input is in degrees (True) or radians (False)
    :return: degrees longitude
    """
    if deg:
        lat = np.deg2rad(lat)

    m_per_deg_lon = 111132.954 * np.cos(lat)
    return m / m_per_deg_lon


def iterative_algorithm(curve_list, curve_array, lat_0, lon_0, az):
    """

    :param curve: array
    :param lat_0:
    :param lon_0:
    :param az: PERPENDICULAR
    :return:
    """

    lat = lat_0
    lon = lon_0

    measured_curve = np.column_stack((curve_array[:, IDX.t], curve_array[:, IDX.f] - curve_array[:, IDX.fb]))
    curve_len = measured_curve.shape[0]

    time_arr = Time(curve_array[:, IDX.t], format="unix")
    r = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
             v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
             v_z=curve_array[:, IDX.vz] * unit.km / unit.s,
             obstime=time_arr)
    r_sat_arr = r.cartesian.without_differentials()
    v_sat_arr = r.velocity

    for i in range(ITER_LIMIT):
        trial_curve = np.empty((curve_len, 2))
        # pick a new location
        lat += m_to_deg_lat(STEP, lat) * np.cos(az)
        lon += m_to_deg_lon(STEP, lat) * np.sin(az)
        alt = 0

        r_user_arr = (EarthLocation.from_geodetic([lon] * curve_len, [lat] * curve_len, [alt] * curve_len)
                      .get_itrs(obstime=time_arr).cartesian.without_differentials())
        r_user_arr: CartesianRepresentation
        r_sat_arr: CartesianRepresentation
        v_sat_arr: CartesianDifferential

        # calculate doppler curve
        rel_vel = v_sat_arr * (r_sat_arr - r_user_arr) / (r_sat_arr - r_user_arr).norm()
        print(rel_vel)



