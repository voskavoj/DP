from src.radio.iridium_frame_operations import decompose_generic_frame
from src.config.setup import WORKING_DIR

EXP_NAME = "val05"
FROM = 2500
TO = 6000

DATA_PATH = WORKING_DIR + f"Data\\{EXP_NAME}\\"


def fits(time):
    if FROM is not None and TO is not None:
        return FROM <= time <= TO
    elif FROM is not None:
        return FROM <= time
    elif TO is not None:
        return time <= TO


with open(DATA_PATH + "decoded.txt", "r") as file:
    frames = file.readlines()

with open(DATA_PATH + "decoded_full.txt", "w") as file:
    for line in frames:
        file.write(line.strip() + "\n")


with open(DATA_PATH + "decoded.txt", "w") as file:
    for line in frames:
        rel_time, freq = decompose_generic_frame(line)
        rel_time /= 1000  # ms to s
        if fits(rel_time):
            file.write(f"{line.strip()}\n")

