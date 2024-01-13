import numpy as np
from astropy.time import Time, TimeDelta
import astropy.units as unit

from src.radio.iridium_channels import map_sat_id_to_tle_id
from src.satellites.predictions import predict_satellite_positions


SATELLITE_FRAME_COUNT_FILTER = 100


def process_received_frames(frames_array: np.array, start_time: str, satellites: dict):
    """
    Process received frames from radio into navigation data
    Input - numpy array of foats: satellite id | relative time | frequency | base frequency
    Ouptut - list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS)

    Prepares absolute time, extracts unique transmitting satellites, assigns them actual ID, predicts satellite positions and composes navigation data.

    :param frames_array: array of received frames
    :param start_time: start time of recording as ISO UTC string
    :param satellites: dicitonary of satellites as Satellite for ONE constellation
    :return: navigation data as list
    """
    # calculate actual times
    array_len = frames_array.shape[0]
    base_times = Time([start_time] * array_len, format="iso", scale="utc")
    rel_times = TimeDelta(frames_array[:, 1] * unit.ms)
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
        except ValueError:
            continue
        nav_data.append([times[i], frames_array[i, 2], frames_array[i, 3], pos_itrs[j, i]])

    return nav_data
