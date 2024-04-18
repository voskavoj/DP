import pickle
import matplotlib.pyplot as plt

from src.config.setup import *
from src.navigation.data_processing import NavDataArrayIndices as IDX, find_curves
from src.utils.plots import ColorPicker


OFFSET = 0


with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
    saved_nav_data = pickle.load(file)
with open(DATA_PATH + TEST_DATA_FILE, "rb") as file:
    test_nav_data = pickle.load(file)

color_picker = ColorPicker()

plt.figure()
plt.title("Actual vs. Simulated data")
plt.xlabel("Unix time [s]")
plt.ylabel("Frequency [Hz]")
plt.figtext(0.15, .15, f"Offset: {OFFSET/1e3:.0f} kHz")


curves = find_curves(test_nav_data)
for curve_array in curves:
    plt.plot(curve_array[:, IDX.t], curve_array[:, IDX.f], marker=".", ms=1.5, alpha=0.25,
             color=color_picker.pick_color(int(curve_array[0, IDX.sat_id])),
             label=str(int(curve_array[0, IDX.sat_id])) + " (sim)")

plt.gca().set_prop_cycle(None)  # reset plot colors

curves = find_curves(saved_nav_data)
for curve_array in curves:
    plt.scatter(curve_array[:, IDX.t], curve_array[:, IDX.f] - OFFSET, marker=".",
                color=color_picker.pick_color(int(curve_array[0, IDX.sat_id])),
                label=str(int(curve_array[0, IDX.sat_id])) + " (act)")

plt.legend()
plt.show()
