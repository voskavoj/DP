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
