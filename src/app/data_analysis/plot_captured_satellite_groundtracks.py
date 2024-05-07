import pickle
import matplotlib.pyplot as plt
from astropy.coordinates import ITRS
from astropy import units as unit
from mpl_toolkits.basemap import Basemap
import numpy as np

from src.config.setup import *
from src.navigation.data_processing import NavDataArrayIndices as IDX, find_curves
from src.utils.plots import ColorPicker
from src.config.locations import LOCATIONS


with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
    saved_nav_data = pickle.load(file)

color_picker = ColorPicker()


curves = find_curves(saved_nav_data, min_curve_length=10, max_time_gap=60)

sat_tracks = list()
min_lon, max_lat, max_lon, min_lat = False, False, False, False
avg_lon, avg_lat = 0, 0
for curve_array in curves:
    print("#", end="")
    sat_track = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
                     v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
                     v_z=curve_array[:, IDX.vz] * unit.km / unit.s).earth_location.geodetic
    min_lon = min(sat_track.lon.min().value, min_lon) if min_lon is not False else sat_track.lon.min().value
    max_lat = max(sat_track.lat.max().value, max_lat) if max_lat is not False else sat_track.lat.max().value
    max_lon = max(sat_track.lon.max().value, max_lon) if max_lon is not False else sat_track.lon.max().value
    min_lat = min(sat_track.lat.min().value, min_lat) if min_lat is not False else sat_track.lat.min().value
    sat_tracks.append((str(int(curve_array[0, IDX.sat_id])), sat_track))
    avg_lon += sat_track.lon.mean().value
    avg_lat += sat_track.lat.mean().value
avg_lat /= len(curves)
avg_lon /= len(curves)
print()
print(min_lat, min_lon, max_lat, max_lon)
print(avg_lat, avg_lon)
home_lon, home_lat = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1]


plt.figure(figsize=(5, 6))
m = Basemap(llcrnrlon=min_lon - 2, llcrnrlat=min_lat - 2, urcrnrlon=max_lon + 2, urcrnrlat=max_lat + 2,  # 6 for legend
            rsphere=(6378137.00, 6356752.3142), resolution='l', projection='merc')

for sat_id, sat_track in sat_tracks:
    print("#", end="")
    m.plot(sat_track.lon, sat_track.lat, ".", ms=1.5, latlon=True,
           # color=color_picker.pick_color(sat_id), label=sat_id)
           color="red", label="_nolegend_")
m.plot([], [], color="red", label="Satellite track")
m.plot(home_lon, home_lat, "o", color="blue", latlon=True, label="User position")
m.plot(avg_lon, avg_lat, "o", color="green", latlon=True, label="Average position")
m.drawcoastlines()
m.fillcontinents()
m.drawcountries()
m.drawrivers(color="blue")
m.drawparallels(np.arange(0, 80, 10), labels=[0, 1, 0, 0])
m.drawmeridians(np.arange(-90, 90, 10), labels=[0, 0, 0, 1])

plt.legend()
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1), ncol=6)
plt.show()
