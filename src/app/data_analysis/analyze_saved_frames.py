import numpy as np
import matplotlib.pyplot as plt

from src.config.setup import *
from src.radio.iridium_offline_radio import IridiumOfflineRadio


GENERATE_NEW_DATA = True

SAVE_FILE = "saved_output_full_decode_drop_2"

# open file
if GENERATE_NEW_DATA:
    with open(DATA_PATH + FRAME_FILE, "r") as file:
        frames = file.readlines()
    radio = IridiumOfflineRadio(frames, file_is_parsed=True)
    frames_array = np.array(radio.get_frames())
    del radio
    np.save(DATA_PATH + SAVE_FILE, frames_array)
else:
    frames_array = np.load(DATA_PATH + SAVE_FILE + ".npy", allow_pickle=True)


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
