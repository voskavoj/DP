import numpy as np
import pickle
from matplotlib import pyplot as plt

from src.radio.iridium_offline_radio import IridiumOfflineRadio
from src.satellites.download_tle import download_tles
from src.navigation.data_processing import process_received_frames
from src.navigation.curve_fit_method import solve

USE_SAVED_DATA = True
SAVED_DATA_FILE = "test_nav_data.pickle"
CONSTELLATIONS = ("Iridium", )
DATA_PATH = "Data\\exp03\\"
FRAME_FILE = "decoded.txt"
START_TIME = "2024-01-07 18:18:29"  # UTC

# ---------------------------- init
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
start_time = START_TIME

if USE_SAVED_DATA:
    with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
        nav_data = pickle.load(file)
else:
    # load frames from radio
    radio = IridiumOfflineRadio(DATA_PATH + FRAME_FILE, file_is_parsed=True)
    # frames_array: satellite ID | relative time | received frequency | base frequency
    frames_array = np.array(radio.get_frames())
    # process nav_data
    # list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS) | ID
    nav_data = process_received_frames(frames_array, start_time, satellites["Iridium"])
    with open(DATA_PATH + SAVED_DATA_FILE, "wb") as file:
        pickle.dump(nav_data, file)

# nav_data = nav_data[1700:2000]
# dopp_data = [(nav_data[i][4], nav_data[i][0].value.split()[1], nav_data[i][2], nav_data[i][2] - nav_data[i][1]) for i in range(len(nav_data))]
# for d in dopp_data:
#     print(*d)
print("Init done.")

# ---------------------------- navigation
solve(nav_data, satellites["Iridium"])
plt.show()

