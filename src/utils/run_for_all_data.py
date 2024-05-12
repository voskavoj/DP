import pickle

from src.config.parameters import CurveFitMethodParameters
from src.config.setup import *
from src.satellites.download_tle import download_tles
from src.utils.data import dump_data, load_data


def run_for_all_data(cb_solve, test_name):
    results = dict()
    for exp_name in VALIDATION_DATA_SETS:
        print(f"Running for {exp_name}")

        data_path = WORKING_DIR + f"Data\\{exp_name}\\"

        with open(data_path + SAVED_DATA_FILE, "rb") as file:
            saved_nav_data = pickle.load(file)

        satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=data_path)
        default_parameters = CurveFitMethodParameters()

        # run callback
        retval = cb_solve(saved_nav_data, satellites["Iridium"], default_parameters, exp_name)
        results[exp_name] = retval

    dump_data("validation\\" + f"{test_name}", results)


def load_results(test_name, index=-1, unpack=False):
    results = load_data("validation\\" + f"{test_name}", index)
    results: dict

    if unpack:
        return list(results.values())
    else:
        return results
