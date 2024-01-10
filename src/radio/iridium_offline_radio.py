from src.radio.iridium_identification import find_doppler_params, identify_transmitting_iridium_satellites_from_ira_frame
from src.radio.iridium_frame_operations import decompose_iridium_frame, is_valid_iridium_frame, get_frame_confidence
from time import monotonic
from astropy.time import Time, TimeDelta
from astropy import units as unit
import numpy as np


"""
    Buffer of frames looks like this:
    satellite ID | time | received frequency | base frequency
"""


CLOSEST_REL_DOPPLER_SHIFT_THRESHOLD = 0.5  # 10%
TRANSMISSION_SET_TIMEOUT = 30  # seconds
MINIMUM_IRA_FRAME_CONFIDENCE = 65  # percent


class IridiumOfflineRadio:
    def __init__(self, file_path: str, start_time: Time, simulate=False):
        with open(file_path, "r") as file:
            self.frames = file.readlines()
        self.start_time = start_time
        self.container = RXFramesContainer()
        self._index = 0

    def get_frames(self, num=100):
        """

        :param num: -1 for all
        :return:
        """
        if self._index >= len(self.frames):
            return False

        if num == -1:
            num = len(self.frames)

        for frame in self.frames[self._index : min(self._index + num, len(self.frames))]:
            if not is_valid_iridium_frame(frame):
                continue

            rel_time, _ = decompose_iridium_frame(frame)
            self.container.clear_old(rel_time / 1000)
            self.container.add_frame(frame)

        self._index += num

        print(f"Time {round(rel_time / 1000)}s | "
              f"Unique TX sets: {len(self.container.content)} | "
              f"Identified sats {[x.sat_id for x in self.container.content if x.sat_id]}")

        return self.container.collect()


class RXFramesContainer:
    def __init__(self):
        self.content: list[TransmissionSet] = list()

    def add_frame(self, frame):
        rel_time, freq = decompose_iridium_frame(frame)
        doppler, rel_doppler, base_freq, channel_id = find_doppler_params(freq)
        if not doppler:
            return

        if not self.content:
            transmission_set = TransmissionSet(rel_doppler)
            self.content.append(transmission_set)
        else:
            closest_tx_set = min(self.content, key=lambda x: abs(x.rel_dopp_shift - rel_doppler))
            closest_rel_doppler_shift = closest_tx_set.rel_dopp_shift

            if abs(rel_doppler - closest_rel_doppler_shift) / closest_rel_doppler_shift < CLOSEST_REL_DOPPLER_SHIFT_THRESHOLD:
                transmission_set = closest_tx_set
                transmission_set.rel_dopp_shift = rel_doppler
            else:
                # print(f"New transmission set with rel_doppler {rel_doppler} vs {closest_rel_doppler_shift}")
                transmission_set = TransmissionSet(rel_doppler)
                self.content.append(transmission_set)

        transmission_set.add_frame(rel_time, freq, base_freq)

        if channel_id == "RA" and not transmission_set.sat_id and get_frame_confidence(frame) >= MINIMUM_IRA_FRAME_CONFIDENCE:
            sat_id, initial_condition = identify_transmitting_iridium_satellites_from_ira_frame(frame)
            if sat_id:
                transmission_set.set_sat_id(sat_id)
                self.add_initial_condition(initial_condition)

    def collect(self):
        out = list()
        for tx_set in self.content:
            out.extend(tx_set.out_buffer)
            tx_set.out_buffer.clear()
        return out

    def clear_old(self, time=None):
        if not time:
            time = monotonic()

        for tx_set in self.content:
            if not tx_set.is_up_to_date(time):
                self.content.remove(tx_set)

    def add_initial_condition(self, initial_condition):
        # todo
        pass


class TransmissionSet:
    def __init__(self, rel_dopp_shift):
        self.sat_id = None
        self.rel_dopp_shift = rel_dopp_shift
        self.last_added_time = 0

        self.no_id_buffer = list()
        self.out_buffer = list()

    def add_frame(self, rel_time, freq, base_freq):
        self.last_added_time = rel_time / 1000
        if self.sat_id:
            self.out_buffer.append(np.array([self.sat_id, rel_time, freq, base_freq], dtype=float))
        else:
            self.no_id_buffer.append([rel_time, freq, base_freq])

    def set_sat_id(self, sat_id):
        self.sat_id = float(sat_id)

        for tx in self.no_id_buffer:
            self.add_frame(*tx)
        self.no_id_buffer.clear()

    def is_up_to_date(self, time=None):
        if not time:
            time = monotonic()
        y = time - self.last_added_time < TRANSMISSION_SET_TIMEOUT
        return time - self.last_added_time < TRANSMISSION_SET_TIMEOUT
