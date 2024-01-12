import subprocess

IRIDIUM_PARSER_PATH = "External/iridium-toolkit/iridium-parser.py"


def parse_raw_iridium_frames(raw_frames: str | list, harder=True) -> list[str]:
    if isinstance(raw_frames, list):
        raw_frames = "\n".join(raw_frames)

    prc = subprocess.run(["python", IRIDIUM_PARSER_PATH, "--harder" if harder else ""],
                         input=raw_frames, text=True, capture_output=True)

    if prc.returncode != 0:
        print(prc.stderr)
        raise ValueError("Error in parsing iridium frames")

    return prc.stdout.split("\n")


def get_frame_confidence(frame: str) -> int:
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


def decompose_ibc_frame(frame: str) -> tuple[int, float, int]:
    """
    Decompose an IBC frame into its components.

    Example parsed frame:
    IBC: p-349207-e000 000044922.3846 1624230912 93% -84.44|-137.62|20.52 138 DL bc:0 sat:057 cell:24 0 slot:1 ...

    :param frame: ibc frame
    :return: sat_id, rel_time, freq
    """
    frame = frame.split()

    rel_time = float(frame[2])
    freq = int(frame[3])

    sat_id = False
    for part in frame[8:]:
        try:
            if "sat" in part:
                sat_id = int(part.split("sat:")[1])
                break
        except ValueError:
            print(f"Part {part} of frame could not be parsed")

    return sat_id, rel_time, freq


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
