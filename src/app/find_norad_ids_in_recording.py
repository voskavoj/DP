"""
    Identify NORAD IDs within Iridium frames
"""

import numpy as np
from astropy.time import Time, TimeDelta
import astropy.units as unit
from sgp4.api import SatrecArray
from collections import defaultdict

from src.config.setup import *
from src.radio.iridium_frame_operations import decompose_ira_frame
from src.satellites.download_tle import download_tles, unpack
from src.satellites.predictions import find_closest_tle_id_to_ira_id_by_position
from src.utils.data import dump_data, load_data
from src.radio.iridium_channels import TLE_ID_TO_IRI_ID_MAP


LOG_LINES = False
LOAD_DATA_IDX = None


def find_norad_ids_in_recording(satellites, frame_file=FRAME_FILE, print_log=False):
    """
        Find the NORAD IDs in the recording
    :param satellites: TLEs
    :param frame_file: decoded frames
    :param print_log: whether to print details
    :return: list of most likely NORAD ID per frame
    """
    sat_list = unpack(satellites)

    with open(DATA_PATH + frame_file, "r") as file:
        frames = file.readlines()

    found_ids = dict()
    ress = list()
    for i, frame in enumerate(frames):
        if not print_log and i % 500 == 0:
            print(f"\r{100 * i / len(frames):.1f}%", end="")
        if frame.startswith("IRA"):
            iri_id, rel_time, freq, ira_data = decompose_ira_frame(frame)
            beam_id, iri_lat, iri_lon, iri_alt, x, y, z = ira_data

            if not (iri_id and rel_time and iri_lat is not False and iri_lon is not False and iri_alt is not False
                    and 700 < iri_alt < 900):
                continue
        else:
            continue

        if not found_ids.get(iri_id):
            found_ids[iri_id] = list()

        # try to identify
        time = Time(START_TIME, format="iso", scale="utc") + TimeDelta(rel_time * unit.ms)
        satrecs = SatrecArray([sat.satrec for sat in sat_list])

        # predict
        tle_id, dist, tle_lat, tle_lon, tle_alt = find_closest_tle_id_to_ira_id_by_position(
            time, satrecs, sat_list, iri_lat, iri_lon, iri_alt)

        if print_log:
            print(f"ID {iri_id:03d} is Iridium {tle_id:03d}, distance {dist / 1000:.0f} km | "
                  f"Frame coordinates: {iri_lat:.2f} deg, {iri_lon:.2f} deg, {iri_alt} km | "
                  f"Pred. coordinates: {tle_lat:.2f}, {tle_lon:.2f}, {tle_alt:.0f}")

        ress.append([iri_id, tle_id, dist, iri_lat, iri_lon, iri_alt, tle_lat, tle_lon, tle_alt])

    return ress


def print_found_norad_ids_in_recording(ress):
    """
        Visualise the found ID map
    :param ress: results of find_norad_ids_in_recording
    """
    ress_arr = np.array(ress)

    # Dictionary to store results
    result = defaultdict(lambda: defaultdict(list))

    # Iterate over rows of ress_arr
    for row in ress_arr:
        iri_id, closest_sat_id, dist = row[0], row[1], row[2]
        result[closest_sat_id][iri_id].append(dist)

    # Sort the result by occurrences
    sorted_result = {}
    for sat_id, iri_data in result.items():
        sorted_iri_data = sorted(iri_data.items(), key=lambda x: len(x[1]), reverse=True)
        sorted_result[sat_id] = sorted_iri_data

    # Print the sorted result
    for sat_id, iri_data in sorted(sorted_result.items()):
        print(f"Iridium {sat_id:.0f}:")
        for iri_id, dist_list in iri_data:
            occurrences = len(dist_list)
            avg_dist = np.mean(dist_list)

            if sat_id not in TLE_ID_TO_IRI_ID_MAP:
                status = "NEWS"
            elif iri_id == TLE_ID_TO_IRI_ID_MAP[sat_id]:
                status = "SAME"
            else:
                status = "DIFF"

            if occurrences > 10 and avg_dist < 100e3:
                dict_style = f" | {sat_id:.0f}: {iri_id:.0f},  # {occurrences}"
            else:
                dict_style = ""
            print(
                f"  IRI ID: {iri_id: 4.0f} ({status}), Occurrences: {occurrences: 4d}, Average Distance: {avg_dist / 1000: 7.1f} km" + dict_style)


if __name__ == "__main__":
    if LOAD_DATA_IDX is None:
        ress = find_norad_ids_in_recording(download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH),
                                           FRAME_FILE, print_log=LOG_LINES)
        dump_data(f"iridium_id_comparison_results_{EXP_NAME}", ress)
    else:
        ress = load_data(f"iridium_id_comparison_results_{EXP_NAME}_{LOAD_DATA_IDX}")

    print_found_norad_ids_in_recording(ress)
