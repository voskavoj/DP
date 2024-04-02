import numpy as np
import pickle

from src.config.setup import *
from src.radio.iridium_offline_radio import IridiumOfflineRadio
from src.satellites.download_tle import download_tles
from src.navigation.data_processing import process_received_frames


# ---------------------------- init
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
start_time = START_TIME

# load frames from radio
radio = IridiumOfflineRadio(DATA_PATH + FRAME_FILE, file_is_parsed=True, non_ira_frames=False, drop_frequencies=True)
# frames_array: satellite ID | relative time | received frequency | base frequency
frames_array = np.array(radio.get_frames())
# process nav_data
# list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS) | ID
nav_data = process_received_frames(frames_array, start_time, satellites["Iridium"])
with open(DATA_PATH + SAVED_DATA_FILE, "wb") as file:
    pickle.dump(nav_data, file)

print("Processing done.")
