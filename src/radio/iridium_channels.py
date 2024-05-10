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
    100: 73,   # 111
    102: 112,  # 24 + 598
    103: 103,  # 103 + 51
    104: 110,  # 53 + 830
    106: 114,  # 147
    107: 115,  # 259
    108: 2,    # 442
    109: 4,    # 23
    110: 9,    # 6 + 315
    111: 16,   # 17 + 508
    112: 17,   # 35 + 769
    113: 24,  # 89
    114: 26,   # 75 + 291
    116: 28,  # 20
    117: 29,  # 77
    118: 33,  # 102
    119: 36,   # 567 + 566
    120: 38,  # 114
    121: 42,  # 152
    122: 44,   # 553 + 416
    123: 48,  # 129
    125: 69,   # 70
    126: 71,  # 139
    128: 78,   # 340
    129: 79,   # 156
    130: 85,  # 86
    131: 87,  # 65
    132: 88,   # 245
    133: 89,   # 27
    134: 92,  # 59
    135: 93,  # 46
    136: 99,   # 775
    137: 104,  # 39
    138: 109,  # 100
    139: 57,   # 735
    141: 51,   # 47
    147: 7,    # 245
    151: 111,  # 150
    152: 22,   # 206
    154: 94,   # 316
    155: 25,   # 422
    156: 46,   # 304 + 544
    158: 18,   # 285 + 540
    159: 49,   # 158 + 464
    160: 90,   # 226 + 555
    163: 3,    # 129 + 389
    164: 13,   # 398
    165: 23,   # 90 + 390
    166: 96,   # 73 + 331
    167: 67,  # 160
    168: 68,  # 67
    171: 81,  # 194
    172: 72,  # 89
    173: 65,  # 58
    180: 50,  # 95

    # 169: 57,   # 275
    # 179: 78,   # 669
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
