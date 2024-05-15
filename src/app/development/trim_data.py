import os

from src.radio.iridium_frame_operations import decompose_generic_frame
from src.config.setup import *

from src.app.process_offline_data import process_offline_data
from src.app.data_analysis.find_offset_and_drift_in_data import find_offset_and_drift_in_data


def fits(time, fr=None, to=None):
    if fr is not None and to is not None:
        return fr <= time <= to
    elif fr is not None:
        return fr <= time
    elif to is not None:
        return time <= to
    else:
        return True


if __name__ == "__main__":
    if not os.path.exists(DATA_PATH + SAVED_DATA_FILE):
        process_offline_data()
    find_offset_and_drift_in_data()

    while True:
        fr = int(input("From time: "))
        to = int(input("To time: "))
        fr = fr if fr != 0 else None
        to = to if to != 0 else None
        print(f"From: {fr}, To: {to}")
        ok = input("OK? (y/n) ")
        if ok == "y":
            break

    with open(DATA_PATH + "decoded.txt", "r") as file:
        frames = file.readlines()

    if not os.path.exists(DATA_PATH + "decoded_full.txt"):
        print("Creating a copy")
        with open(DATA_PATH + "decoded_full.txt", "w") as file:
            for line in frames:
                file.write(line.strip() + "\n")
    else:
        print("Full copy already exists")

    print("Trimming data")
    with open(DATA_PATH + "decoded.txt", "w") as file:
        for line in frames:
            rel_time, freq = decompose_generic_frame(line)
            rel_time /= 1000  # ms to s
            if fits(rel_time, fr, to):
                file.write(f"{line.strip()}\n")

    print("Processing trimmed data")
    process_offline_data()
