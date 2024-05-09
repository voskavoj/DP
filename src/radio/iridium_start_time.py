import numpy as np
from astropy.time import Time

from src.radio.iridium_frame_operations import decompose_ibc_frame, get_frame_confidence
from src.config.setup import *


def compute_start_time(frames, confidence_threshold=80, time_spread_threshold=10):
    """
    Compute start time of the recording based on the received IBC frames

    :param frames: list of decoded frames (as text)
    :param confidence_threshold: minimum confidence level in %
    :param time_spread_threshold: maximum time spread in ms/min
    :return: start time (unix), correction factor (s/s)
    """
    print("Parsing frames")
    time_list = list()
    for i, frame in enumerate(frames):
        if i % 100 == 0:
            print(f"\r{100 * i / len(frames):.1f}%", end="")
        if frame.startswith("IBC"):
            _, rel_time, _, sync_time = decompose_ibc_frame(frame)
            if sync_time is False:
                continue

            conf = get_frame_confidence(frame)
            if conf < confidence_threshold:
                continue

            rel_time /= 1000  # to seconds

            start_time = sync_time - rel_time
            time_list.append([rel_time, sync_time, start_time])

            # print(f"Rel time: {rel_time:.0f} | Start time: {Time(start_time, format='unix', scale='utc').iso}")
    print(f"\rDone")

    # convert to numpy array
    time_arr = np.array(time_list)
    rel_time_arr = time_arr[:, 0]
    start_time_arr = time_arr[:, 2]
    duration = rel_time_arr[-1] - rel_time_arr[0]

    # filter out outliers
    start_time_med = np.median(start_time_arr)
    outlier_mask = np.abs(start_time_arr - start_time_med) < time_spread_threshold

    rel_time_arr = rel_time_arr[outlier_mask]
    start_time_arr = start_time_arr[outlier_mask]

    # compute start time and correction factor
    if start_time_arr.shape[0] > 100:
        start_time = np.median(start_time_arr[:100])
        start_time_end = np.median(start_time_arr[-100:])
        corr_factor = 1 + (start_time_end - start_time) / duration
    else:
        start_time = np.median(start_time_arr)
        start_time_end = None
        corr_factor = 1

    # plot
    # plt.figure()
    # p = np.poly1d(fit)
    # ptime = np.linspace(np.min(rel_time_arr), np.max(rel_time_arr), 2)
    # plt.plot(rel_time_arr, start_time_arr)
    # plt.plot(ptime, p(ptime), color="black", label="Fit")
    # plt.figtext(.15, .12, f"Fit parameters ($s/s^n$): {fit}")
    # plt.show()

    print(f"Start time goes from {Time(start_time, format="unix", scale="utc").iso} "
          f"to {Time(start_time_end, format="unix", scale="utc").iso} in {duration:.0f} s ({duration/60:.0f} min) - "
          f"correction factor is 1 + {(corr_factor - 1)*1e6:.3f} us/s")

    return start_time, corr_factor


if __name__ == "__main__":
    with open(DATA_PATH + FRAME_FILE, "r") as file:
        frames = file.readlines()
    start_time, corr_factor = compute_start_time(frames)
    with open(DATA_PATH + "start_time.txt", "w") as file:
        frames = file.write(f"{Time(start_time, format="unix", scale="utc").iso}\n")
