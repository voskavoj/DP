import numpy as np
import matplotlib.pyplot as plt

from src.radio.iridium_offline_radio import IridiumOfflineRadio


GENERATE_NEW_DATA = True

DATA_PATH = "Data\\exp03\\"
# FRAME_FILE = "output.bits"
FRAME_FILE = "decoded.txt"
SAVED_DATA_FILE = "saved_output_full_decode_drop_2"

# open file
if GENERATE_NEW_DATA:
    radio = IridiumOfflineRadio(DATA_PATH + FRAME_FILE, file_is_parsed=True)
    frames_array = None

    while frames := radio.get_frames():
        if frames_array is None:
            frames_array = np.array(frames)
        else:
            frames_array = np.append(frames_array, frames, axis=0)

    del radio
    np.save(DATA_PATH + SAVED_DATA_FILE, frames_array)
else:
    frames_array = np.load(DATA_PATH + SAVED_DATA_FILE + ".npy", allow_pickle=True)


print(frames_array.shape)
x = frames_array[:, 1] / 1000
y = frames_array[:, 2]
b = frames_array[:, 3]
c = frames_array[:, 0]

# plt.figure()
# plt.plot(frames_array[:, 2], ".")
# color = frames_array[:, 0]
# rgb = plt.get_cmap('jet')(color)
# plt.plot(frames_array[:, 2], ".", color=rgb)
plt.scatter(x, y, s=1, c=c, cmap='tab10')
plt.scatter(x, b, s=0.2, c=c, cmap='tab20')


# plt.plot(frames_array[:, 3], ".")
#
# plt.figure()
# plt.plot(frames_array[:, 0], ".")
plt.show()
