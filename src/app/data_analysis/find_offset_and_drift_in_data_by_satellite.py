import pickle

import numpy as np
from astropy.coordinates import ITRS
import astropy.units as unit
from matplotlib import pyplot as plt

from src.config.locations import LOCATIONS
from src.config.setup import *
from src.navigation.data_processing import find_curves, NavDataArrayIndices as IDX
from src.satellites.download_tle import download_tles
from src.navigation.curve_fit_method import solve, check_trial_curve

satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
lon, lat, alt = LOCATIONS[LOCATION][0], LOCATIONS[LOCATION][1], LOCATIONS[LOCATION][2]


with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
    nav_data = pickle.load(file)
detected_curves = find_curves(nav_data)

results = list()
for curve_array in detected_curves:
    r = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
             v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
             v_z=curve_array[:, IDX.vz] * unit.km / unit.s)
    measured_curve = np.column_stack(
        (curve_array[:, IDX.t], curve_array[:, IDX.f] - curve_array[:, IDX.fb], curve_array[:, IDX.fb]))
    r_sat_arr = np.column_stack((curve_array[:, IDX.x], curve_array[:, IDX.y], curve_array[:, IDX.z]))
    v_sat_arr = np.column_stack((curve_array[:, IDX.vx], curve_array[:, IDX.vy], curve_array[:, IDX.vz]))

    off_range = range(-20000, +20000, 100)
    res = list()
    for off in off_range:
        sos = check_trial_curve(lat, lon, alt, off, 0, measured_curve, r_sat_arr, v_sat_arr)
        res.append((sos, off))
    off = min(res, key=lambda x: x[0])[1]

    off_range = range(off-100, off+100, 1)
    res = list()
    for off in off_range:
        sos = check_trial_curve(lat, lon, alt, off, 0, measured_curve, r_sat_arr, v_sat_arr)
        res.append((sos, off))
    off = min(res, key=lambda x: x[0])[1]

    sat_id = int(curve_array[0, IDX.sat_id])
    time_start = curve_array[0, IDX.t]
    time_end = curve_array[-1, IDX.t]
    results.append((time_start, time_end, off, sat_id))
    print(f"Offset for {curve_array[0, IDX.sat_id]:.0f}: {off} Hz")

plt.figure()
results = sorted(results, key=lambda x: x[0])
min_time = min(results, key=lambda x: x[0])[0]
for time_start, time_end, off, sat_id in results:
    t = time_start - min_time + (time_end - time_start) / 2
    plt.scatter(t, off, label=str(sat_id))

data = np.array(results)
data[:, 0] -= min_time
data[:, 1] -= min_time
data[:, 0] = data[:, 0] + (data[:, 1] - data[:, 0]) / 2
fit = np.polyfit(data[:, 0], data[:, 2], 1)
p = np.poly1d(fit)
print("Fit: ", fit)

plt.plot(data[:, 0], p(data[:, 0]), color="black", label="Fit")
plt.figtext(.15, .12, f"Fit parameters ($Hz/s^n$): {fit}")
plt.legend()
plt.show()

