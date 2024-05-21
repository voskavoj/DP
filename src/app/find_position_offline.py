"""
    Find position from processed data
    This is the main development script

    Currently set to run the curve fit method
"""

import pickle

import numpy as np
from astropy.coordinates import EarthLocation
from matplotlib import pyplot as plt

from src.config.locations import LOCATIONS
from src.config.setup import *
from src.navigation.calculations import latlon_distance
from src.satellites.download_tle import download_tles
from src.navigation.curve_fit_method import solve, C
from src.config.parameters import CurveFitMethodParameters


LON_HOME, LAT_HOME, ALT_HOME = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]


def slice_data_into_time_chunks(data, chunk_time):
    if chunk_time is None:
        return [data]

    time_chunks = list()
    start_time, end_time = data[0, 0], data[-1, 0]
    i = 0

    while start_time + chunk_time * i < end_time:
        mask = (data[:, 0] >= start_time + chunk_time * i) & (data[:, 0] < start_time + chunk_time * (i + 1))
        time_chunks.append(data[mask])
        i += 1

    return time_chunks


with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
    saved_nav_data = pickle.load(file)

# with open(DATA_PATH + TEST_DATA_FILE, "rb") as file:
#     test_nav_data = pickle.load(file)

satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)["Iridium"]
parameters = CurveFitMethodParameters()


est_state = None
for data in slice_data_into_time_chunks(saved_nav_data, None):
    res = solve(data, satellites, parameters, init_state=est_state)
    lat, lon, alt, off, dft = res
    # est_state = res
    print(f"Estimated state: dist {latlon_distance(LAT_HOME, lat, LON_HOME, lon):04.0f} m | lat {lat:.3f}, lon {lon:.3f}, alt {alt:.0f}, off {off:.0f}, dft {dft:.3f}")

plt.show()
