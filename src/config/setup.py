# Experiment settings
EXP_NAME = "val01"
LOCATION = "HOME"
START_TIME = None
VALIDATION_DATA_SETS = ["val01", "val02", "val03", "val04", "val05"]

# Paths
WORKING_DIR = "C:\\Git\\Personal\\DP\\"
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

    dump_results = True


if START_TIME is None:
    with open(DATA_PATH + "start_time.txt", "r") as file:
        START_TIME = file.readlines()[0].strip()
