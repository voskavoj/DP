import subprocess


IRIDIUM_PARSER_PATH = "GNURadio/iridium-toolkit/iridium-parser.py"


def is_valid_iridium_frame(frame: str) -> bool:
    return frame.startswith("RAW:") and len(frame) >= 150


def decompose_iridium_frame(frame: str) -> tuple[float, int]:
    """
    Decompose a (non-parsed) iridium frame into its components.

    Example frame:
    RAW: i-349207-t1 0000922.8643 1626149248 N:33.10-87.82 I:00000000000 100% 0.06553 432 00110...
    0    | 1                            | 2              | 3         | 4                | 5        | 6          | 7            | 8      | 9
    type | unix timestamp of rec. start | ms since start | frequency | signal magnitude | frame ID | confidence | signal level | length | demodulated bits

    :param frame: non-parsed iridium frame from gr-iridium as a string
    :return:
    """
    frame = frame.split()

    rel_time = float(frame[2])
    freq = int(frame[3])

    return rel_time, freq


def get_frame_confidence(frame: str) -> int:
    return int(frame.split()[6].replace("%", ""))


def get_unix_timestamp_from_iridium_frame(frame: str) -> int:
    """
    Get unix timestamp from an iridium frame.

    :param frame: non-parsed iridium frame from gr-iridium as a string
    :return: unix timestamp as int
    """
    return int(frame.split()[1].split("-")[1])


def parse_and_decompose_ira_frames(frames: list | str) -> list[tuple[int, int, float, float, int, int, int, int] | None]:
    """
    Parse iridium frames from a string or a list.

    Example parsed frame:
    IRA: p-289693-e000 000025108.8066 1626314368  91% -29.47|-090.82|30.79 107 DL sat:067 beam:27 xyz=(+1390,+0092,+1121) pos=(+38.82/+003.79) alt=797 RAI:48 ?00 bc_sb:28 P00: {TRUNCATED}

    :param frames: iridium frames
    :return: parsed frames
    """

    if isinstance(frames, list):
        frames = "\n".join(frames)

    prc = subprocess.run(["python", IRIDIUM_PARSER_PATH, "--harder"],
                         input=frames, text=True, capture_output=True)

    if prc.returncode != 0:
        print("Error in parsing iridium frames")
        print(prc.stdout)
        print(prc.stderr)
        return None

    # filter IRA frames
    parsed_ira_frames = [line for line in prc.stdout.split("\n") if line.startswith("IRA")]
    decomposed_ira_frames = list()

    if not parsed_ira_frames:
        print("No IRA frames found")
    else:
        for ira_frame in parsed_ira_frames:
            if decomposed := decompose_parsed_ira_frame(ira_frame):
                decomposed_ira_frames.append(decomposed)

    return decomposed_ira_frames


def decompose_parsed_ira_frame(ira_frame: str) -> tuple[int, int, float, float, int, int, int, int]:
    """
    Decompose an ira frame into its components.

    Example parsed frame:
    IRA: p-289693-e000 000025108.8066 1626314368  91% -29.47|-090.82|30.79 107 DL sat:067 beam:27 xyz=(+1390,+0092,+1121) pos=(+38.82/+003.79) alt=797 RAI:48 ?00 bc_sb:28 P00: {TRUNCATED}

    :param ira_frame: ira frame
    :return: sat_id, beam_id, lat, lon, alt, x, y, z; any of which can be False
    """
    ira_frame = ira_frame.split()
    sat_id, beam_id, lat, lon, alt, x, y, z = [False] * 8

    for part in ira_frame:
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

    return sat_id, beam_id, lat, lon, alt, x, y, z
