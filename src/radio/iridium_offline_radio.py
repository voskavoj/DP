import copy
import subprocess

from src.radio.iridium_identification import find_tx_base_frequency
from src.radio.iridium_frame_operations import decompose_parsed_ira_frame
from astropy.time import Time
import numpy as np

"""
    Buffer of frames looks like this:
    satellite ID | time | received frequency | base frequency
"""

IRIDIUM_PARSER_PATH = "GNURadio/iridium-toolkit/iridium-parser.py"

CLOSEST_REL_DOPPLER_SHIFT_THRESHOLD = 0.5  # 10%
TRANSMISSION_SET_TIMEOUT = 30  # seconds
MINIMUM_IRA_FRAME_CONFIDENCE = 65  # percent


class IridiumOfflineRadio:
    def __init__(self, file_path: str, file_is_parsed=False):
        with open(file_path, "r") as file:
            frames = file.read()

        self.frames_list = list()

        if not file_is_parsed:
            print("Parsing frames")
            prc = subprocess.run(["python", IRIDIUM_PARSER_PATH],
                                 input=frames, text=True, capture_output=True)
            frames = prc.stdout

        frames = frames.split("\n")

        print("Analyzing.")
        for frame in frames:
            if frame.startswith("IRA") or frame.startswith("IBC"):
                sat_id, beam_id, lat, lon, alt, x, y, z = decompose_parsed_ira_frame(frame)
                split_frame = frame.split()
                rel_time = float(split_frame[2])
                freq = int(split_frame[3])
                base_freq = find_tx_base_frequency(freq, drop=False)

                if sat_id is not False and base_freq is not False:
                    self.frames_list.append(np.array([sat_id, rel_time, freq, base_freq]))

    def get_frames(self, num=100):
        """

        :param num: -1 for all
        :return:
        """
        if self.frames_list:
            retval = copy.deepcopy(self.frames_list)
            self.frames_list.clear()
            return retval
        else:
            return False
