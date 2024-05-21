import numpy as np
import pickle

from astropy.time import Time

from src.config.setup import *
from src.radio.iridium_offline_radio import IridiumOfflineRadio
from src.radio.iridium_start_time import compute_start_time
from src.satellites.download_tle import download_tles
from src.navigation.data_processing import process_received_frames
from src.utils.data import save_data


def process_offline_data(start_time=None, save_file=SAVED_DATA_FILE, save=True):
    """
    Process the demodulated and decoded frames into navigation data

    #) Load TLEs from Data directory (arg. offline_dir)
    #) Load demodulated frames from Data directory
    #) Find start time by the use IBC frame - radio/iridium_start_time
    #) Process frames into array of time, f, fb, sat_id - radio/iridium_offline_radio
    #) Process frames into navigation data - navigation/data_processing

    :param start_time: time of start of the recording
    :param save_file: name of the file to save the data to
    :param save: whether to save the data
    :return: processed navigation data
    """
    satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
    with open(DATA_PATH + FRAME_FILE, "r") as file:
        frames = file.readlines()

    # find start time
    if start_time is None:
        start_time, time_corr_factor = compute_start_time(frames)
        start_time = Time(start_time, format="unix", scale="utc").iso
    else:
        time_corr_factor = 1

    # load frames from radio
    radio = IridiumOfflineRadio(frames, file_is_parsed=True, non_ira_frames=False, drop_frequencies=True)
    # frames_array: satellite ID | relative time | received frequency | base frequency
    frames_array = np.array(radio.get_frames())
    # process nav_data
    # list: absolute time (Time) | frequency (float) | base frequency (float) | satellite position at time (ITRS) | ID
    nav_data = process_received_frames(frames_array, start_time, satellites["Iridium"],
                                       time_correction_factor=time_corr_factor)

    if save:
        with open(DATA_PATH + save_file, "wb") as file:
            pickle.dump(nav_data, file)

    print("Processing done.")
    return nav_data


if __name__ == "__main__":
    process_offline_data()
