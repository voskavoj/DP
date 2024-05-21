"""
    Generate test data for the navigation algorithm based on the satellites within the measured data
    Parametrised below
"""

import numpy as np
import pickle
from astropy.time import Time, TimeDelta

from src.config.setup import *
from src.config.locations import LOCATIONS
from src.radio.iridium_channels import map_tle_id_to_sat_id, IRA_BASE_FREQUENCY
from src.satellites.download_tle import download_tles
from src.navigation.data_processing import process_received_frames, find_curves, NavDataArrayIndices as IDX
from src.satellites.predictions import predict_satellite_doppler_shift

TEST_SAT_MIN_SAMPLES = 50
BASE_FREQ = IRA_BASE_FREQUENCY
TEST_DATA_MINUTES = 160
SKIPPED_MINUTES = 100
DECIMATION = 4
DOPP_LIMIT = 36000

# ---------------------------- init
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
start_time = Time(START_TIME) + TimeDelta(SKIPPED_MINUTES * 60, format="sec")
lon, lat, alt = LOCATIONS[LOCATION]

# ---------------------------- get satellites present in actual data
with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
    saved_nav_data = pickle.load(file)


test_sat_list = list()
for curve_array in find_curves(saved_nav_data):
    sat = satellites["Iridium"][str(int(curve_array[0, IDX.sat_id]))]
    sat_no, length = sat.number, curve_array.shape[0]
    if length > TEST_SAT_MIN_SAMPLES:
        test_sat_list.append(sat)
        print(f"\t {sat.name}: {length: 5d} samples")

# [rel_time, shifted frequency, abs_doppler_shift, rel_doppler_shift]
dopp_predictions = predict_satellite_doppler_shift(test_sat_list, lat, lon, alt, BASE_FREQ,
                                                   start_time, duration=TEST_DATA_MINUTES * 60, step=0.1)

# [sat_id, rel_time (ms), freq, base_freq]
frames_list = list()
for i in range(len(test_sat_list)):
    sat_id = map_tle_id_to_sat_id(test_sat_list[i].number)
    if sat_id is False:
        raise ValueError(i, test_sat_list[i], sat_id)
    data_len = dopp_predictions.shape[1]

    sat_id_arr = np.ones((data_len, 1)) * float(sat_id)
    base_freq_arr = np.ones((data_len, 1)) * float(BASE_FREQ)
    dopp_preds = np.hstack((sat_id_arr, dopp_predictions[i, :, 0:1] * 1000, dopp_predictions[i, :, 1:2], base_freq_arr))

    for j, pred in enumerate(dopp_preds):
        if DECIMATION and j % DECIMATION != 0:
            continue
        if pred[1] == 0 or pred[2] == 0:
            continue
        frames_list.append(pred)

# filter empty frames
# for i, fr in enumerate(frames_list):
#     if fr[1] == 0 or fr[2] == 0:
#         frames_list.pop(i)

# satellite ID | relative time | received frequency | base frequency
frames_array = np.array(frames_list)
print(frames_array.shape)

nav_data = process_received_frames(frames_array, START_TIME, satellites["Iridium"], time_correction_factor=1)
print(len(nav_data))

with open(DATA_PATH + TEST_DATA_FILE, "wb") as file:
    pickle.dump(nav_data, file)
