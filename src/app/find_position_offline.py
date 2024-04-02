import pickle
from matplotlib import pyplot as plt

from src.config.setup import *
from src.satellites.download_tle import download_tles
from src.navigation.curve_fit_method import solve


def find_position_from_offline_data(satellites, data_filename):
    with open(DATA_PATH + data_filename + ".pickle", "rb") as file:
        nav_data = pickle.load(file)

    solve(nav_data, satellites["Iridium"])


satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
find_position_from_offline_data(satellites, "test_nav_data")
find_position_from_offline_data(satellites, "saved_nav_data")

plt.show()
