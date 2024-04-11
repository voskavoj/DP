import pickle

import numpy as np
from matplotlib import pyplot as plt

from src.config.setup import *
from src.navigation.calculations import add_gaussian_noise_and_offset
from src.satellites.download_tle import download_tles
from src.navigation.curve_fit_method import solve
from src.navigation.data_processing import NavDataArrayIndices as IDX


# test variations - mean, std, offset
VARIATIONS = [
    (0, 0, 0),
    (0, 0, -15000),
    (0, 100, 0),
    (0, 0, 20000),
]

satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
with open(DATA_PATH + "test_nav_data" + ".pickle", "rb") as file:
    nav_data = pickle.load(file)

for mean, std, offset in VARIATIONS:
    print(f"Mean: {mean}, Std: {std}, Offset: {offset}")

    data = np.copy(nav_data)
    data[:, IDX.f] = add_gaussian_noise_and_offset(nav_data[:, IDX.f],
                                                   noise_mean=mean, noise_std_dev=std, offset=offset)
    solve(data, satellites["Iridium"])


plt.show()
