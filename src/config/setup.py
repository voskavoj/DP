# Experiment settings
EXP_NAME = "exp04"
LOCATION = "HOME"
START_TIME = None

# Paths
WORKING_DIR = "C:\\Git\\Personal\\DP\\"
TMP_PATH = WORKING_DIR + "tmp\\"
DATA_PATH = WORKING_DIR + f"Data\\{EXP_NAME}\\"

# Global files, ...
CONSTELLATIONS = ("Iridium",)
FRAME_FILE = "decoded.txt"
SAVED_DATA_FILE = "saved_nav_data.pickle"
TEST_DATA_FILE = "test_nav_data.pickle"

if START_TIME is None:
    with open(DATA_PATH + "start_time.txt", "r") as file:
        START_TIME = file.readlines()[0].strip()
