import os
import pickle

from src.config.setup import TMP_PATH

FILE_IDX = 0
FIG_IDX = 0


def dump_data(name, data, overwrite=False):
    global FILE_IDX

    create_tmp_dir()
    filename = f"{TMP_PATH}{name}_{FILE_IDX}.pickle"

    while True:
        if not overwrite and os.path.exists(filename):
            FILE_IDX += 1
            filename = f"{TMP_PATH}{name}_{FILE_IDX}.pickle"
        else:
            break

    print(f"Saving data as {name}_{FILE_IDX}")
    with open(filename, "wb") as file:
        pickle.dump(data, file)


def load_data(name):
    filename = f"{TMP_PATH}{name}.pickle"
    with open(filename, "rb") as file:
        data = pickle.load(file)
    return data


def get_fig_filename(name):
    global FIG_IDX

    create_tmp_dir()
    filename = f"{TMP_PATH}{name}_{FIG_IDX}.png"

    while True:
        if os.path.exists(filename):
            FIG_IDX += 1
            filename = f"{TMP_PATH}{name}_{FIG_IDX}.png"
        else:
            return filename


def create_tmp_dir():
    if os.path.exists(TMP_PATH):
        return
    else:
        os.makedirs(TMP_PATH)


# def clear_tmp_dir():
#     for file in os.listdir(TMP_PATH):
#         os.remove(TMP_PATH + file)
