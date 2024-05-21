import sys
import subprocess
from matplotlib import pyplot as plt
import numpy as np
from astropy.time import Time

from src.radio.iridium_offline_radio import IridiumOfflineRadio
from src.radio.iridium_start_time import compute_start_time
from src.navigation.data_processing import process_received_frames
from src.satellites.download_tle import download_tles
from src.navigation.curve_fit_method import solve
from src.config.parameters import CurveFitMethodParameters


IRIDIUM_PARSER_PATH = "External/iridium-toolkit/iridium-parser.py"


def main(data_path, tle_path=None):
    """
    Main function to find the position from offline data

    :param data_path: path to the demodulated frames
    :param tle_path: path to the TLE file
    """
    # decode frames
    print("Decoding frames, this may take a while")
    print(f"Command: python {IRIDIUM_PARSER_PATH} -p {data_path} --harder")
    prc = subprocess.run(["python", IRIDIUM_PARSER_PATH, "-p", data_path, "--harder"],
                         text=True, capture_output=True)
    if prc.returncode != 0:
        print(prc.stderr)
        raise ValueError("Error in parsing iridium frames")
    frames = prc.stdout.split("\n")

    # download TLEs
    if not tle_path:
        satellites = download_tles(constellations=("Iridium",))
    else:
        print("Using provided TLES")
        satellites = download_tles(constellations=("Iridium",), offline_dir=tle_path)

    # find start time
    start_time, time_corr_factor = compute_start_time(frames)
    start_time = Time(start_time, format="unix", scale="utc").iso

    # process frames
    radio = IridiumOfflineRadio(frames, file_is_parsed=True, non_ira_frames=False, drop_frequencies=True)
    frames_array = np.array(radio.get_frames())
    nav_data = process_received_frames(frames_array, start_time, satellites["Iridium"],
                                       time_correction_factor=time_corr_factor)

    # solve user position
    parameters = CurveFitMethodParameters()
    res = solve(nav_data, satellites["Iridium"], parameters)
    lat, lon, alt, off, dft = res
    print(f"Estimated state: lat {lat:.3f}, lon {lon:.3f}, alt {alt:.0f}, off {off:.0f}, dft {dft:.3f}")

    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python find_position_from_offline_data.py <path to demodulated frames> <directory with tle file (optional)>")
        print("Example: python find_position_from_offline_data.py /Data/exp01/output.bits" "/Data/exp01/")
        print(f"iridium-toolkit needs to be downloaded at {IRIDIUM_PARSER_PATH}")
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("Too many arguments")
