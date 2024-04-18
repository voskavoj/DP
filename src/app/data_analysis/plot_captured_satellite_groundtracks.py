import pickle
import matplotlib.pyplot as plt
from astropy.coordinates import ITRS
from astropy import units as unit
from mpl_toolkits.basemap import Basemap
import numpy as np

from src.config.setup import *
from src.navigation.data_processing import NavDataArrayIndices as IDX, find_curves
from src.utils.plots import ColorPicker
from src.navigation.calculations import latlon_distance
from src.navigation.data_processing import nav_data_to_array
from src.utils.data import load_data, get_fig_filename
from src.config.locations import LOCATIONS


with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
    saved_nav_data = pickle.load(file)

color_picker = ColorPicker()


curves = find_curves(saved_nav_data, min_curve_length=10, max_time_gap=60)

sat_tracks = list()
min_lon, max_lat, max_lon, min_lat = False, False, False, False
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
print()
print(min_lat, min_lon, max_lat, max_lon)

home_lon, home_lat = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1]


plt.figure()
m = Basemap(llcrnrlon=min_lon - 2, llcrnrlat=min_lat - 2, urcrnrlon=max_lon + 2, urcrnrlat=max_lat + 2,
            rsphere=(6378137.00, 6356752.3142), resolution='l', projection='merc')

for sat_id, sat_track in sat_tracks:
    print("#", end="")
    m.plot(sat_track.lon, sat_track.lat, ".", ms=1.5, latlon=True,
           color=color_picker.pick_color(sat_id), label=sat_id)
m.drawcoastlines()
m.fillcontinents()
m.drawcountries()
m.drawrivers(color="blue")
m.drawparallels(np.arange(40, 65, 5), labels=[0, 1, 0, 0])
m.drawmeridians(np.arange(0, 30, 5), labels=[0, 0, 0, 1])


plt.legend()
plt.show()
