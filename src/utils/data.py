import os
import pickle

from src.config.setup import TMP_PATH, DATA_PATH

FIG_IDX = 0


def dump_data(name, data):
    """
        Dump data to tmp folder

    :param name: filename, suffixed by index
    :param data: data
    """
    create_tmp_dir()

    file_idx = 0
    filename = f"{TMP_PATH}{name}_{file_idx}.pickle"

    while os.path.exists(filename):
        file_idx += 1
        filename = f"{TMP_PATH}{name}_{file_idx}.pickle"

    print(f"Saving data as {filename}")
    with open(filename, "wb") as file:
        pickle.dump(data, file)


def save_data(name, data):
    """
        Save data to Data folder, overwrites existing file
    :param name: filename
    :param data: data
    """

    filename = f"{DATA_PATH}{name}.pickle"
    with open(filename, "wb") as file:
        pickle.dump(data, file)


def load_data(name, index=None):
    """
    Load data from a dump.

    :param name: filename
    :param index: file index - if None, none is inserted, if -1, the last existing index is used
    :return: data from dump
    """
    if index is None:
        filename = f"{TMP_PATH}{name}.pickle"
    elif index == -1:
        index = 0
        while os.path.exists(f"{TMP_PATH}{name}_{index+1}.pickle"):
            index += 1
        filename = f"{TMP_PATH}{name}_{index}.pickle"
    elif isinstance(index, int):
        filename = f"{TMP_PATH}{name}_{index}.pickle"
    else:
        raise ValueError("Index must be an integer or None")

    try:
        with open(filename, "rb") as file:
            data = pickle.load(file)
    except FileNotFoundError:
        print(f"File {filename} not found")
        return None
    return data


def get_fig_filename(name, idx=True):
    """
    Get next available figure filename
    :param name: name
    :param idx: use index, if not, overriden
    :return: filename
    """
    global FIG_IDX

    create_tmp_dir()
    if not idx:
        return f"{TMP_PATH}{name}.png"

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
