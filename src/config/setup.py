"""
    Setup of code

    This does NOT contain algorithm parameters, which are in config/parameters

    Important path is the WORKING_DIR, should point to the root of the project

    EXP_NAME is the name of the experiment, which is used to load data from the Data folder
    LOCATION is the location of the experiment, defined in locations, used to calculate distance to est. position
    START_TIME is the start time of the experiment, if None, it is loaded from the data folder
    VALIDATION_DATA_SETS is a list of all the validation data sets

    DEBUG is a class containing all the debug flags, used to control the output of the program
"""

raise NotImplementedError("EXP_NAME, LOCATION and WORKING_DIR need to be set in src/config/setup.py. Remove this line afterward.")

# Experiment settings
EXP_NAME = ""
LOCATION = ""
START_TIME = None
VALIDATION_DATA_SETS = ["val01", "val02", "val03", "val04", "val05", "val06", "val07", "val08", "val09", "val10"]

# Paths
WORKING_DIR = ""
TMP_PATH = WORKING_DIR + "tmp\\"
DATA_PATH = WORKING_DIR + f"Data\\{EXP_NAME}\\"

# Global files, ...
CONSTELLATIONS = ("Iridium",)
FRAME_FILE = "decoded.txt"
SAVED_DATA_FILE = "saved_nav_data.pickle"
TEST_DATA_FILE = "test_nav_data.pickle"


# LOGGING
class DEBUG:
    log_detected_curves = False
    log_detail_progress = False

    plot_results = False
    plot_analyzed_curve = False
    plot_final_curve_fit = False

    dump_results = False


if START_TIME is None:
    try:
        with open(DATA_PATH + "start_time.txt", "r") as file:
            START_TIME = file.readlines()[0].strip()
    except FileNotFoundError:
        START_TIME = None
        print("WARNING: Start time not found in experiment data folder!")
