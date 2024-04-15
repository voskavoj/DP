import copy
import numpy as np

from src.radio.iridium_channels import find_tx_base_frequency, IRA_BASE_FREQUENCY
from src.radio.iridium_frame_operations import decompose_ira_frame, decompose_ibc_frame, decompose_generic_frame, \
    parse_raw_iridium_frames

"""
    Buffer of frames looks like this:
    satellite ID | relative time | received frequency | base frequency
"""


class IridiumOfflineRadio:
    def __init__(self, file_path: str, file_is_parsed=False, non_ira_frames=True, drop_frequencies=True):
        with open(file_path, "r") as file:
            frames = file.readlines()

        self.frames_list = list()

        if not file_is_parsed:
            frames = parse_raw_iridium_frames(frames)

        print("Analyzing.")
        for frame in frames:
            if frame.startswith("IRA"):
                sat_id, rel_time, freq, ira_data = decompose_ira_frame(frame)
                base_freq = IRA_BASE_FREQUENCY
            elif non_ira_frames and frame.startswith("IBC"):
                sat_id, rel_time, freq = decompose_ibc_frame(frame)
                base_freq = find_tx_base_frequency(freq, drop=drop_frequencies)
            else:
                continue

            if sat_id and base_freq is not False:
                self.frames_list.append(np.array([sat_id, rel_time, freq, base_freq]))

    def get_frames(self):
        return copy.deepcopy(self.frames_list)
