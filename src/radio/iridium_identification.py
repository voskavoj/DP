import subprocess

from src.radio.iridium_frame_operations import parse_and_decompose_ira_frames
from src.satellites.satellite import generate_sattelite_id

IRIDIUM_PARSER_PATH = "GNURadio/iridium-toolkit/iridium-parser.py"
IRA_BASE_FREQUENCY = 1626270800
_CHANNEL_SPACING = (1 / 3) * 10e5 / 8
_DUPLEX_FREQS = [round(1616000000 + _CHANNEL_SPACING / 2 + _CHANNEL_SPACING * i) for i in range(240)]
_SIMPLEX_FREQ = [1626104200, 1626145800, 1626270800, 1626395800, 1626437500]
CHANNEL_FREQS = _DUPLEX_FREQS + _SIMPLEX_FREQ

CHANNEL_IDS = {ch: v for v, ch in enumerate(_DUPLEX_FREQS)}
_SIMPLEX_IDS = ["P1", "P2", "RA", "P3", "P4"]
# _SIMPLEX_IDS = [300, 310, 320, 330, 340]
CHANNEL_IDS.update({ch: v for ch, v in zip(_SIMPLEX_FREQ, _SIMPLEX_IDS)})
CHANNEL_ID_THRESHOLD = _CHANNEL_SPACING


def find_doppler_params(freq):
    """
    Find doppler shift, relative doppler shift, base frequency and channel ID from a frequency.

    :param freq: frequency (int)
    :return: doppler_shift, relative_doppler_shift, base_frequency, channel_id
    """
    base_freq = find_tx_base_frequency(freq, drop=True)
    if not base_freq:
        return False, False, False, False
    channel_id = CHANNEL_IDS[base_freq]
    doppler = freq - base_freq
    rel_doppler = (freq - base_freq) / base_freq

    return doppler, rel_doppler, base_freq, channel_id


def find_tx_base_frequency(freq, drop=False):
    base_freq = min(CHANNEL_FREQS, key=lambda x: abs(x - freq))
    if drop and abs(base_freq - freq) > CHANNEL_ID_THRESHOLD:
        return False
        # print(f"WARNING: Frequency {freq} is outside of channel threshold for channel {CHANNEL_IDS[base_freq]} "
        #       f"by {round((abs(base_freq - freq) - CHANNEL_ID_THRESHOLD) * 1e-3, 1)} kHz")
    return base_freq


def find_doppler_shift(freq):
    base_freq = find_tx_base_frequency(freq)
    return freq - base_freq


def find_channel_id(freq):
    return CHANNEL_IDS[find_tx_base_frequency(freq, drop=False)]

#
# def identify_transmitting_iridium_satellites_from_ira_frame(ira_frame, require_one=False):
#     decomposed_ira_frames = _parse_ira_frames(frames)
#     if not decomposed_ira_frames:
#         return None
#
#     tx_sats = dict()
#     for x in decomposed_ira_frames:
#         tx_sats[x[0]] = x
#     tx_sats = list(tx_sats.values())
#
#     if require_one and len(tx_sats) != 1:
#         print(f"Found {len(tx_sats)} transmitting satellites")
#         return False
#     elif require_one:
#         return tx_sats[0]
#     else:
#         return tx_sats


def identify_transmitting_iridium_satellites_from_ira_frame(ira_frame):
    decomposed_ira_frames = parse_and_decompose_ira_frames(ira_frame)
    if not decomposed_ira_frames:
        return False, False
    else:
        sat_id, beam_id, lat, lon, alt, x, y, z = decomposed_ira_frames[0]
        if sat_id == 65:
            print(ira_frame)
        return get_iridium_id(sat_id), (lat, lon, alt, x, y, z)


def get_iridium_id(sat_id: int):
    # todo identification based on inner number
    return generate_sattelite_id("Iridium", sat_id)

