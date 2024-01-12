import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, TEME, ITRS
from astropy.time import Time, TimeDelta
import astropy.units as unit
from astropy.visualization import time_support
from sgp4.api import SatrecArray

from src.radio.iridium_channels import map_sat_id_to_tle_id
from src.radio.iridium_offline_radio import IridiumOfflineRadio
from src.satellites.download_tle import download_tles

CONSTELLATIONS = ("Iridium", )
DATA_PATH = "Data\\exp03\\"
FRAME_FILE = "decoded.txt"
START_TIME = "2024-01-07 18:17:40"  # UTC

# ---------------------------- init
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
radio = IridiumOfflineRadio(DATA_PATH + FRAME_FILE, file_is_parsed=True)

start_time = START_TIME

# ---------------------------- loop (here once)
# load frames from radio
# frame = [satellite_id, time, freq, base_freq]
frames_array = np.array(radio.get_frames())

# predict satellite positions at each time


# get actual times
array_len = frames_array.shape[0]
base_times = Time([start_time] * array_len, format="iso", scale="utc")
rel_times = TimeDelta(frames_array[:, 1] * unit.ms)
times = base_times + rel_times

# filter unique satellites
tx_satellites, count = np.unique(frames_array[:, 0], return_counts=True)
tx_satellites = tx_satellites[count > 100]  # filter out satellites with less than 100 frames

# get actual satellite IDs
tx_satellites_ids = [str(map_sat_id_to_tle_id(sat_id)) for sat_id in tx_satellites]
tx_satrecs = [satellites["Iridium"][sat_id].satrec for sat_id in tx_satellites_ids]

# SGP4
tx_satrec_array = SatrecArray(tx_satrecs)
jds = np.array(times.jd1)
jfs = np.array(times.jd2)

errcodes, positions, velocities = tx_satrec_array.sgp4(jds, jfs)
pos_teme = CartesianRepresentation(x=positions[:, :, 0], y=positions[:, :, 1], z=positions[:, :, 2], unit=unit.km)
vel_teme = CartesianDifferential(velocities[:, :, 0], velocities[:, :, 1], velocities[:, :, 2], unit=unit.km / unit.s)
r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=times)
pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))


# nav array: time, freq, base_freq, sat_pos_itrs
nav_data_list = list()

for i in range(array_len):
    try:
        j = tx_satellites_ids.index(str(map_sat_id_to_tle_id(frames_array[i, 0])))
    except ValueError:
        continue
    nav_data_list.append([times[i], frames_array[i, 2], frames_array[i, 3], pos_itrs[j, i]])

print("Done.")







