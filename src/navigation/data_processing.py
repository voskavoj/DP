import numpy as np
from astropy.time import Time, TimeDelta
import astropy.units as unit

from src.radio.iridium_channels import map_sat_id_to_tle_id
from src.satellites.predictions import predict_satellite_positions


SATELLITE_FRAME_COUNT_FILTER = 100
TIME_CORRECTION_FACTOR = 1  # ms/ms


class NavDataIndices:
    abs_time = 0
    t = abs_time
    freq = 1
    f = freq
    base_freq = 2
    fb = base_freq
    sat_pos = 3
    pos = sat_pos
    sat_id = 4
    id = sat_id


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

    for i in range(array_len):
        try:
            j = tx_satellites_ids.index(str(map_sat_id_to_tle_id(frames_array[i, 0])))
            sat = map_sat_id_to_tle_id(frames_array[i, 0])
        except ValueError:
            continue
        nav_data.append([times[i], frames_array[i, 2], frames_array[i, 3], pos_itrs[j, i], sat])

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
                             sat_pos.x.value, sat_pos.y.value, sat_pos.z.value,
                             sat_pos.v_x.value, sat_pos.v_y.value, sat_pos.v_z.value]

    return nav_data_array


def find_curves(nav_data_array):
    IDX = NavDataArrayIndices
    # sort nav_data_array by sat_id, base_freq, and abs_time
    sorted_indices = np.lexsort((nav_data_array[:, IDX.t], nav_data_array[:, IDX.fb], nav_data_array[:, IDX.id]))

    sorted_nav_data = nav_data_array[sorted_indices]

    # find unique combinations of sat_id and base_freq
    unique_combinations, indices = np.unique(sorted_nav_data[:, [IDX.id, IDX.fb]], axis=0, return_index=True)

    split_arrays = []

    # Split sorted_nav_data into arrays based on unique combinations
    for i in range(len(indices)):
        if i == len(indices) - 1:
            split_array = sorted_nav_data[indices[i]:]
        else:
            split_array = sorted_nav_data[indices[i]:indices[i + 1]]

        # verify the splitting works
        assert np.all(split_array[:, IDX.fb] == split_array[0, IDX.fb])  # is one base frequency
        assert np.all(split_array[:, IDX.id] == split_array[0, IDX.id])  # is one satellite ID
        assert np.all(np.diff(split_array[:, IDX.t]) >= 0)  # is sorted by time

        if True:  # todo is_valid_curve(split_array):
            split_arrays.append(split_array)

    return split_arrays
