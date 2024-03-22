IRIDIUM_PARSER_PATH = "GNURadio/iridium-toolkit/iridium-parser.py"
IRA_BASE_FREQUENCY = 1626270800
_CHANNEL_SPACING = (1 / 3) * 10e5 / 8
_DUPLEX_FREQS = [round(1616000000 + _CHANNEL_SPACING / 2 + _CHANNEL_SPACING * i) for i in range(240)]
_SIMPLEX_FREQ = [1626104200, 1626145800, 1626270800, 1626395800, 1626437500]
CHANNEL_FREQS = _DUPLEX_FREQS + _SIMPLEX_FREQ

CHANNEL_IDS = {ch: v for v, ch in enumerate(_DUPLEX_FREQS)}
_SIMPLEX_IDS = ["P1", "P2", "RA", "P3", "P4"]
CHANNEL_IDS.update({ch: v for ch, v in zip(_SIMPLEX_FREQ, _SIMPLEX_IDS)})
CHANNEL_ID_THRESHOLD = _CHANNEL_SPACING


def find_tx_base_frequency(freq, drop=False):
    base_freq = min(CHANNEL_FREQS, key=lambda x: abs(x - freq))
    if drop and abs(base_freq - freq) > CHANNEL_ID_THRESHOLD:
        return False
        # print(f"WARNING: Frequency {freq} is outside of channel threshold for channel {CHANNEL_IDS[base_freq]} "
        #       f"by {round((abs(base_freq - freq) - CHANNEL_ID_THRESHOLD) * 1e-3, 1)} kHz")
    return base_freq


def find_channel_id(freq):
    return CHANNEL_IDS[find_tx_base_frequency(freq, drop=False)]


TLE_ID_TO_IRI_ID_MAP = {
    102: 112,  # 24
    103: 103,  # 103
    104: 110,  # 53
    109: 4,    # 23
    110: 9,    # 6
    111: 16,   # 17
    112: 17,   # 35
    114: 26,   # 75
    119: 36,   # 567
    122: 44,   # 553
    156: 46,   # 304
    158: 18,   # 285
    159: 49,   # 158
    160: 90,   # 226
    163: 3,    # 129
    165: 23,   # 90
    166: 96,   # 73
    169: 57,   # 275
    179: 78,   # 669
}

IRI_ID_TO_TLE_ID_MAP = {iri_id: tle_id for tle_id, iri_id in TLE_ID_TO_IRI_ID_MAP.items()}


if len(IRI_ID_TO_TLE_ID_MAP) != len(TLE_ID_TO_IRI_ID_MAP):
    raise ValueError(f"Duplicate values in IRI_ID_TO_TLE_ID_MAP: "
                     f"{set([x for x in TLE_ID_TO_IRI_ID_MAP.values() if list(TLE_ID_TO_IRI_ID_MAP.values()).count(x) > 1])}")


def map_sat_id_to_tle_id(sat_id: str | int | float) -> int | bool:
    """
    Maps the satellite ID from an IRA/IBC frame to the TLE ID

    :param sat_id: satellite ID from IRA/IBC frame as number
    :return: Satellite ID from TLE as integer
    """
    return IRI_ID_TO_TLE_ID_MAP[int(sat_id)]


def map_tle_id_to_sat_id(tle_id: str | int | float) -> int | bool:
    """
    Maps the satellite ID from a TLE to the ID from IRA/IBC frame

    :param tle_id: Satellite ID from TLE as number
    :return: Satellite ID from IRA/IBC frame as integer
    """
    return TLE_ID_TO_IRI_ID_MAP[int(tle_id)]
