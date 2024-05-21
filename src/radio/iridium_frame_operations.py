"""
    Iriridium frame operations
"""

import subprocess

from astropy.time import Time, TimeDelta
import astropy.units as unit

IRIDIUM_PARSER_PATH = "External/iridium-toolkit/iridium-parser.py"


def parse_raw_iridium_frames(raw_frames: str | list, harder=True) -> list[str]:
    """
        Parse raw iridium frames using the iridium-parser script.
        Make sure the path is correct

    :param raw_frames: list of raw iridium frames or a single frame
    :param harder: use --harder flag
    :return: list of decoded iridium frames
    """
    if isinstance(raw_frames, list):
        raw_frames = "\n".join(raw_frames)

    prc = subprocess.run(["python", IRIDIUM_PARSER_PATH, "--harder" if harder else ""],
                         input=raw_frames, text=True, capture_output=True)

    if prc.returncode != 0:
        print(prc.stderr)
        raise ValueError("Error in parsing iridium frames")

    return prc.stdout.split("\n")


def get_frame_confidence(frame: str) -> int:
    """
    Get the confidence of a decoded iridium frame.
    :param frame: decoded iridium frame
    :return: confidence as int in %
    """
    return int(frame.split()[4].replace("%", ""))


def get_unix_timestamp_from_iridium_frame(frame: str) -> int:
    """
    Get unix timestamp from an iridium frame.

    :param frame: non-parsed iridium frame from gr-iridium as a string
    :return: unix timestamp as int
    """
    return int(frame.split()[1].split("-")[1])


def decompose_ira_frame(frame: str) -> tuple[int, float, int, tuple[int, float, float, int, int, int, int]]:
    """
    Decompose an IRA frame into its components.

    Example parsed frame:
    IRA: p-289693-e000 000025108.8066 1626314368  91% -29.47|-090.82|30.79 107 DL sat:067 beam:27 xyz=(+1390,+0092,+1121) pos=(+38.82/+003.79) alt=797 RAI:48 ?00 bc_sb:28 P00: {TRUNCATED}

    :param frame: ira frame
    :return: sat_id, rel_time, freq, ira_data (beam_id, lat, lon, alt, x, y, z; any of which can be False)
    """
    frame = frame.split()
    
    rel_time = float(frame[2])
    freq = int(frame[3])
    
    sat_id, beam_id, lat, lon, alt, x, y, z = [False] * 8

    for part in frame[8:]:
        try:
            if "sat" in part:
                sat_id = int(part.split("sat:")[1])
            elif "beam" in part:
                beam_id = int(part.split(":")[1])
            elif "xyz" in part:
                x, y, z = part.replace("xyz=(", "").replace(")", "").split(",")
                x, y, z = int(x), int(y), int(z)
            elif "pos" in part:
                lat, lon = part.replace("pos=(", "").replace(")", "").split("/")
                lat, lon = float(lat), float(lon)
            elif "alt" in part:
                alt = int(part.split("alt=")[1])
        except ValueError:
            print(f"Part {part} of frame could not be parsed")
        
    ira_data = (beam_id, lat, lon, alt, x, y, z)

    return sat_id, rel_time, freq, ira_data


def decompose_ibc_frame(frame: str) -> tuple[int, float, int, float]:
    """
    Decompose an IBC frame into its components.

    Example parsed frame:
    IBC: p-349207-e000 000044922.3846 1624230912 93% -84.44|-137.62|20.52 138 DL bc:0 sat:057 cell:24 0 slot:1
        sv_blkn:0 aq_cl:1111111111111111 aq_sb:28 aq_ch:2 00 0000 time:2024-03-17T15:24:12.89Z ...

    :param frame: ibc frame
    :return: sat_id, rel_time (ms), freq (Hz), sync_time (unix)
    """
    frame = frame.split()

    rel_time = float(frame[2])
    freq = int(frame[3])

    sat_id, sync_time, slot = False, False, False
    for part in frame[8:]:
        try:
            if "sat" in part:
                sat_id = int(part.split("sat:")[1])
            if "slot" in part:
                slot = int(part.split("slot:")[1])
            if "time" in part:
                sync_time = Time(part.split("time:")[1], format="isot", scale="utc")
                break  # the time comes after the sat_id
        except ValueError:
            print(f"Part {part} of frame could not be parsed")

    if sync_time and slot:
        sync_time = _correct_ibc_time(sync_time, slot).to_value("unix")
    else:
        sync_time = False

    return sat_id, rel_time, freq, sync_time


def _correct_ibc_time(ibc_time: Time, slot: int):
    """
    Correct the time of an IBC frame for slot number and beginning of frame and signal

    Returns time in the same reference frame as the gr-iridium time

    The time in the IBC frame is the time at the start of the 90 ms L-Band frame that contained the IBC frame
    Source: https://github.com/muccc/iridium-toolkit/issues/40
    Code from: https://github.com/muccc/iridium-toolkit/blob/master/iridiumtk/reassembler/ppm.py#L48-L66

    :param ibc_time: timestamp of the frame (Time)
    :param slot: slot of the frame (int)
    :return: corrected timestamp
    """
    # correct for slot:
    # 1st vs. 4th slot is 3 * (downlink + guard)
    ibc_time += TimeDelta((slot * (3 * float(8.28 + 0.1))) * unit.ms)

    # correct to beginning of frame:
    # guard + simplex + guard + 4*(uplink + guard) + extra_guard
    ibc_time += TimeDelta((1 + 20.32 + 1.24 + 4 * float(8.28 + 0.22) + 0.02) * unit.ms)

    # correct to beginning of signal:
    # gr-iridium timestamp is "the middle of the first symbol of the 12-symbol BPSK Iridium sync word"
    # so correct for 64 symbols preamble & one half symbol.
    ibc_time += TimeDelta((64.5 / 25000) * unit.ms)

    return ibc_time


def decompose_generic_frame(frame: str) -> tuple[float, int]:
    """
    Decompose a generic frame into its components.

    Example parsed frame:
    IBC: p-349207-e000 000044922.3846 1624230912 93% -84.44|-137.62|20.52 138 DL bc:0 sat:057 cell:24 0 slot:1 ...

    :param frame: generic frame
    :return: rel_time, freq
    """
    frame = frame.split()

    rel_time = float(frame[2])
    freq = int(frame[3])

    return rel_time, freq
