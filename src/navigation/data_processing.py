import numpy as np
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
import astropy.units as unit

from src.radio.iridium_channels import map_sat_id_to_tle_id
from src.satellites.predictions import predict_satellite_positions
from src.config.locations import LOCATIONS
from src.config.setup import LOCATION

LON_HOME, LAT_HOME, ALT_HOME = LOCATIONS[LOCATION][0], LOCATIONS[LOCATION][1], LOCATIONS[LOCATION][2]
C = 299792458  # m/s

SATELLITE_FRAME_COUNT_FILTER = 10
TIME_CORRECTION_FACTOR = 1  # ms/ms


class NavDataArrayIndices:
    rel_time = 0
    t = rel_time
    freq = 1
    f = freq
    base_freq = 2
    fb = base_freq
    sat_id = 3
    id = sat_id
    x = 4
    y = 5
    z = 6
    vx = 7
    vy = 8
    vz = 9


IDX = NavDataArrayIndices


def process_received_frames(frames_array: np.array, start_time: str, satellites: dict,
                            time_correction_factor: float = TIME_CORRECTION_FACTOR):
    """
    Process received frames from radio into navigation data
    Input - numpy array of foats: satellite id | relative time | frequency | base frequency
    Ouptut - list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS) | ID (TLE)

    Prepares absolute time, extracts unique transmitting satellites, assigns them actual ID, predicts satellite positions and composes navigation data.

    :param frames_array: array of received frames
    :param start_time: start time of recording as ISO UTC string
    :param satellites: dictionary of satellites as Satellite for ONE constellation
    :param time_correction_factor: Factor to correct radio clock
    :return: navigation data as list
    """
    # calculate actual times
    array_len = frames_array.shape[0]
    base_times = Time([start_time] * array_len, format="iso", scale="utc")
    rel_times = TimeDelta(frames_array[:, 1] * time_correction_factor * unit.ms)
    times = base_times + rel_times

    # get all unique transmitting satellites
    tx_satellites, count = np.unique(frames_array[:, 0], return_counts=True)
    tx_satellites = tx_satellites[count > SATELLITE_FRAME_COUNT_FILTER]  # filter out satellites with few frames

    # get actual satellite IDs
    tx_satellites_ids = [str(map_sat_id_to_tle_id(sat_id)) for sat_id in tx_satellites]
    tx_satrecs = [satellites[sat_id].satrec for sat_id in tx_satellites_ids]

    # SGP4
    pos_itrs = predict_satellite_positions(tx_satrecs, times)

    # nav_data: time, freq, base_freq, sat_pos_itrs
    nav_data = list()
    not_found_iri_ids = list()

    for i in range(array_len):
        try:
            j = tx_satellites_ids.index(str(map_sat_id_to_tle_id(frames_array[i, 0])))
            sat = map_sat_id_to_tle_id(frames_array[i, 0])
        except (ValueError, KeyError):
            not_found_iri_ids.append(int(frames_array[i, 0]))
            continue
        nav_data.append([times[i], frames_array[i, 2], frames_array[i, 3], pos_itrs[j, i], sat])

    if not_found_iri_ids:
        not_found_iri_ids = sorted(list(zip(*np.unique(not_found_iri_ids, return_counts=True))), key=lambda x: x[1], reverse=True)
        print(f"Some Iridium IDs were not found in channel map: " + ",".join([f"{iri_id:03d} ({count}x) " for iri_id, count in not_found_iri_ids]))

    return nav_data_to_array(nav_data)  # todo actual array above


def nav_data_to_array(nav_data: list) -> np.array:
    """
    Convert navigation data to numpy array
    :param nav_data: list of navigation data
    :return: numpy array
    """
    nav_data_array = np.empty((len(nav_data), 10))
    for i, data in enumerate(nav_data):
        abs_time, freq, base_freq, sat_pos, sat_id = data
        nav_data_array[i] = [float(abs_time.to_value("unix")), freq, base_freq, sat_id,
                             sat_pos.x.to("km").value, sat_pos.y.to("km").value, sat_pos.z.to("km").value,
                             sat_pos.v_x.to("km/s").value, sat_pos.v_y.to("km/s").value, sat_pos.v_z.to("km/s").value]

    return nav_data_array


def find_curves(nav_data_array, max_time_gap=0, min_curve_length=0):
    """
    Find curves in navigation data array based on max_time_gap and min_curve_length
    One curve is one base freq. and one sat_id

    :param nav_data_array: array of navigation data
    :param max_time_gap: maximum time gap between consecutive data points in one curve
    :param min_curve_length: minimum length of curve in datapoints
    :return: list of arrays, each with a curve
    """
    # sort nav_data_array by sat_id, base_freq, and abs_time
    sorted_indices = np.lexsort((nav_data_array[:, IDX.t], nav_data_array[:, IDX.fb], nav_data_array[:, IDX.id]))

    sorted_nav_data = nav_data_array[sorted_indices]

    # find unique combinations of sat_id and base_freq
    unique_combinations, indices = np.unique(sorted_nav_data[:, [IDX.id, IDX.fb]], axis=0, return_index=True)

    split_arrays = list()

    # split sorted_nav_data into arrays based on unique combinations
    for i in range(len(indices)):
        if i == len(indices) - 1:
            split_array = sorted_nav_data[indices[i]:]
        else:
            split_array = sorted_nav_data[indices[i]:indices[i + 1]]

        # split the curves based on max_time_gap
        if max_time_gap > 0:
            gap_split_curves = list()
            # find indices where consecutive abs_time values exceed max_time_gap
            gap_indices = np.where(np.diff(split_array[:, 0]) > max_time_gap)[0]

            # split the array at split_indices
            if len(gap_indices) > 0:
                prev_idx = 0
                for idx in gap_indices:
                    gap_split_curves.append(split_array[prev_idx:idx + 1])
                    prev_idx = idx + 1
                gap_split_curves.append(split_array[gap_indices[-1] + 1:])
            else:
                gap_split_curves.append(split_array)
        else:
            gap_split_curves = [split_array]

        for curve in gap_split_curves:
            if min_curve_length and len(curve) < min_curve_length:
                continue

            # verify the splitting works
            assert np.all(curve[:, IDX.fb] == curve[0, IDX.fb])  # is one base frequency
            assert np.all(curve[:, IDX.id] == curve[0, IDX.id])  # is one satellite ID
            assert np.all(np.diff(curve[:, IDX.t]) >= 0)  # is sorted by time
            assert np.all(np.diff(curve[:, IDX.t]) <= max_time_gap) if max_time_gap > 0 else True  # max_time_gap
            assert len(curve) >= min_curve_length if min_curve_length > 0 else True  # min_curve_length

            split_arrays.append(curve)

    return split_arrays


def generate_trial_curve(nav_data, lat=LAT_HOME, lon=LON_HOME, alt=ALT_HOME, off=0, dft=0, abs_freq=True):
    """
    Generate a trial curve based on the input navigation data and the provided parameters.
    Copied from curve_fit_method

    :param nav_data: navigation data array
    :param lat: latitude
    :param lon: longitude
    :param alt: altitude
    :param off: offset
    :param dft: drift
    :param abs_freq: generate absolute frequency (if False, Doppler shift curve is generated)
    :return: trial curve as array
    """

    # prepare data
    curve_array = nav_data
    measured_curve = np.column_stack((curve_array[:, IDX.t],
                                      curve_array[:, IDX.f] - curve_array[:, IDX.fb],
                                      curve_array[:, IDX.fb]))
    r_sat_arr = np.column_stack((curve_array[:, IDX.x], curve_array[:, IDX.y], curve_array[:, IDX.z]))
    v_sat_arr = np.column_stack((curve_array[:, IDX.vx], curve_array[:, IDX.vy], curve_array[:, IDX.vz]))

    # generate trial curve
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
    if abs_freq:
        trial_curve[:, 1] = f_d + f_b + off + f_drift
    else:
        trial_curve[:, 1] = f_d + off + f_drift

    return trial_curve

