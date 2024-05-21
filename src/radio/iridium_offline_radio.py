import copy
import numpy as np

from src.radio.iridium_channels import find_tx_base_frequency, IRA_BASE_FREQUENCY, CHANNEL_ID_THRESHOLD
from src.radio.iridium_frame_operations import decompose_ira_frame, decompose_ibc_frame, decompose_generic_frame, \
    parse_raw_iridium_frames

"""
    Buffer of frames looks like this:
    satellite ID | relative time | received frequency | base frequency
"""


class IridiumOfflineRadio:
    """
        Iririum offline radio class
        Works over saved data
        Written to simplify transition to Real Time calc.

        Call get_frames() to get the processed frames
    """
    def __init__(self, frames: list, file_is_parsed=False, non_ira_frames=True, drop_frequencies=True):
        """
        :param frames: frames
        :param file_is_parsed: if frames are decoded (True) or just demodulated (False)
        :param non_ira_frames: Use non-ira frames too
        :param drop_frequencies: drop frames with frequencies that are outside of channel spacing
        """
        self.frames_list = list()

        if not file_is_parsed:
            frames = parse_raw_iridium_frames(frames)

        print("Analyzing.")
        for frame in frames:
            if frame.startswith("IRA"):
                sat_id, rel_time, freq, ira_data = decompose_ira_frame(frame)
                base_freq = IRA_BASE_FREQUENCY
                if drop_frequencies and abs(base_freq - freq) > CHANNEL_ID_THRESHOLD:
                    continue
            elif non_ira_frames and frame.startswith("IBC"):
                sat_id, rel_time, freq, sync_time = decompose_ibc_frame(frame)
                base_freq = find_tx_base_frequency(freq, drop=drop_frequencies)
            else:
                continue

            if sat_id and base_freq is not False:
                self.frames_list.append(np.array([sat_id, rel_time, freq, base_freq]))

    def get_frames(self):
        return copy.deepcopy(self.frames_list)
