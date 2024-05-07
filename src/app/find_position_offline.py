import pickle

import numpy as np
from matplotlib import pyplot as plt

from src.config.setup import *
from src.satellites.download_tle import download_tles
from src.navigation.curve_fit_method import solve
from src.config.parameters import CurveFitMethodParameters


def find_position_from_offline_data(satellites, data_filename, params):
    with open(DATA_PATH + data_filename, "rb") as file:
        nav_data = pickle.load(file)

    solve(nav_data, satellites["Iridium"], params)


satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
parameters = CurveFitMethodParameters()
find_position_from_offline_data(satellites, SAVED_DATA_FILE, parameters)
# find_position_from_offline_data(satellites, TEST_DATA_FILE, parameters)

plt.show()
