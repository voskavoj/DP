import pickle

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time


DATA_PATH = "Data\\exp03\\"
FRAME_FILE = "decoded.txt"
SAVED_DATA_FILE = "saved_nav_data_w_id.pickle"
TEST_DATA_FILE = "test_nav_data.pickle"
BASE_FREQ = 1626270800
START_TIME = "2024-01-07 18:18:29"  # UTC


def convert_nav_data(data, filter_freq=BASE_FREQ):
    start_time = Time(START_TIME)
    reduced_data = list()
    for d in data:
        if filter_freq and d[2] != filter_freq:
            continue
        reduced_data.append([(d[0] - start_time).sec, d[1], d[2], d[4]])
    return np.array(reduced_data)


with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
    saved_nav_data = convert_nav_data(pickle.load(file))
with open(DATA_PATH + TEST_DATA_FILE, "rb") as file:
    test_nav_data = convert_nav_data(pickle.load(file))

print(saved_nav_data.shape, test_nav_data.shape)
sx = saved_nav_data[:, 0]
sy = saved_nav_data[:, 1]
sb = saved_nav_data[:, 2]
sc = saved_nav_data[:, 3]

tx = test_nav_data[:, 0]
ty = test_nav_data[:, 1]
tb = test_nav_data[:, 2]
tc = test_nav_data[:, 3]

bb = np.ones((len(tx),)) * BASE_FREQ

plt.figure()
plt.scatter(tx, bb, s=0.1, c="red")
plt.scatter(sx, sy, s=1, c=sc, cmap='tab20b')
plt.scatter(tx, ty, s=0.1, c=tc, cmap='tab20b')

plt.show()
