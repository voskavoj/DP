import pickle

import numpy as np
from matplotlib import pyplot as plt

from src.navigation.calculations import add_gaussian_noise_and_offset
from src.satellites.download_tle import download_tles
from src.navigation.curve_fit_method import solve
from src.navigation.data_processing import NavDataArrayIndices as IDX


CONSTELLATIONS = ("Iridium", )
DATA_PATH = "Data\\exp04\\"


with open(DATA_PATH + "test_nav_data" + ".pickle", "rb") as file:
    nav_data = pickle.load(file)
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)

solve(nav_data, satellites["Iridium"])

data = np.copy(nav_data)
data[:, IDX.f] = add_gaussian_noise_and_offset(nav_data[:, IDX.f], noise_mean=0, noise_std_dev=10, offset=0)
solve(data, satellites["Iridium"])

data = np.copy(nav_data)
data[:, IDX.f] = add_gaussian_noise_and_offset(nav_data[:, IDX.f], noise_mean=0, noise_std_dev=100, offset=0)
solve(data, satellites["Iridium"])

data = np.copy(nav_data)
data[:, IDX.f] = add_gaussian_noise_and_offset(nav_data[:, IDX.f], noise_mean=0, noise_std_dev=0, offset=20000)
solve(data, satellites["Iridium"])

plt.show()
