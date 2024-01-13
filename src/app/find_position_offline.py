import numpy as np

from src.radio.iridium_offline_radio import IridiumOfflineRadio
from src.satellites.download_tle import download_tles
from src.navigation.data_processing import process_received_frames

CONSTELLATIONS = ("Iridium", )
DATA_PATH = "Data\\exp03\\"
FRAME_FILE = "decoded.txt"
START_TIME = "2024-01-07 18:17:40"  # UTC

# ---------------------------- init
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
radio = IridiumOfflineRadio(DATA_PATH + FRAME_FILE, file_is_parsed=True)

start_time = START_TIME

# ---------------------------- loop (here once)
# load frames from radio
frames_array = np.array(radio.get_frames())
# process nav_data
# list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS)
nav_data = process_received_frames(frames_array, start_time, satellites["Iridium"])

print(nav_data)
print("Done.")







