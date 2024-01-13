import numpy as np
import pickle

from src.radio.iridium_offline_radio import IridiumOfflineRadio
from src.satellites.download_tle import download_tles
from src.navigation.data_processing import process_received_frames

USE_SAVED_DATA = True
SAVED_DATA_FILE = "saved_nav_data.pickle"
CONSTELLATIONS = ("Iridium", )
DATA_PATH = "Data\\exp03\\"
FRAME_FILE = "decoded.txt"
START_TIME = "2024-01-07 18:17:40"  # UTC

# ---------------------------- init
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
start_time = START_TIME

if USE_SAVED_DATA:
    with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
        nav_data = pickle.load(file)
else:
    # load frames from radio
    radio = IridiumOfflineRadio(DATA_PATH + FRAME_FILE, file_is_parsed=True)
    frames_array = np.array(radio.get_frames())
    # process nav_data
    # list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS)
    nav_data = process_received_frames(frames_array, start_time, satellites["Iridium"])
    with open(DATA_PATH + SAVED_DATA_FILE, "wb") as file:
        pickle.dump(nav_data, file)
print("Init done.")

# ---------------------------- navigation


