"""
    App
    Predicts close satellite passes
"""
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.visualization import time_support
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from src.config.locations import LOCATIONS
from src.satellites.download_tle import download_tles, unpack
from src.satellites.predictions import predict_satellite_visibility

time_support()

LOCATION = "HOME"
CONSTELLATIONS = ("Iridium", "Orbcomm", "Globalstar")

PREDICTION_TIME = 12 * 60 * 60
PREDICTION_STEP = 30
PREDICTION_START = "2023-12-17 08:00:00"
SIM_VIS_ELE_LIM = 10
COLORS = ("tab:red", "tab:orange", "tab:green")


observer_position = EarthLocation.from_geodetic(*LOCATIONS[LOCATION])
satellites_dict = download_tles(CONSTELLATIONS)
satellites_list = unpack(satellites_dict)
time = Time(PREDICTION_START)

prediction_lists = predict_satellite_visibility(satellites_list, observer_position, time,
                                                duration=PREDICTION_TIME, step=PREDICTION_STEP,
                                                elevation_limit=None, log=False, native_output=True)


# fig, ax = plt.subplots()
plt.figure()
# elevation in time
for sat, pred in zip(satellites_list, prediction_lists):
    if not pred:
        continue

    times, elevs, azims = zip(*pred)

    if "Iridium" in sat.name:
        color = COLORS[0]
    elif "Orbcomm" in sat.name:
        color = COLORS[1]
    elif "Globalstar" in sat.name:
        color = COLORS[2]
    else:
        raise KeyError(sat.name)

    plt.plot(Time(times), elevs, color=color)

custom_lines = [Line2D([0], [0], color=COLORS[0], lw=1),
                Line2D([0], [0], color=COLORS[1], lw=1),
                Line2D([0], [0], color=COLORS[2], lw=1)]
plt.ylim(0, 100)
plt.grid()
plt.legend(custom_lines, ['Iridium Next', 'Orbcomm', 'Globalstar'])
plt.ylabel("Elevation (°)")

# simultaneous number
simul_vis = dict()
for sat, pred in zip(satellites_list, prediction_lists):
    if not pred:
        continue
    times, elevs, azims = zip(*pred)
    for t in times:
        simul_vis[t] = (0, 0, 0)

for sat, pred in zip(satellites_list, prediction_lists):
    if not pred:
        continue
    times, elevs, azims = zip(*pred)
    for t, e in zip(times, elevs):
        if e > SIM_VIS_ELE_LIM:
            svi, svo, svg = simul_vis[t]

            if "Iridium" in sat.name:
                svi += 1
            elif "Orbcomm" in sat.name:
                svo += 1
            elif "Globalstar" in sat.name:
                svg += 1
            else:
                raise KeyError(sat.name)

            simul_vis[t] = (svi, svo, svg)

times, svs = zip(*sorted(simul_vis.items()))
svis, svos, svgs = zip(*svs)

plt.figure()
plt.grid(axis='y')
plt.stackplot(Time(times), svis, svos, svgs)

plt.legend(['Iridium Next', 'Orbcomm', 'Globalstar'])
plt.ylabel("Number of simultaneously visible satellites")
plt.show()
