import numpy as np
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, ITRS, TEME, EarthLocation
from astropy.time import Time, TimeDelta
import astropy.units as unit
from sgp4.api import SatrecArray

from src.config.setup import *
from src.radio.iridium_frame_operations import decompose_ira_frame
from src.satellites.download_tle import download_tles, unpack


# ---------------------------- init
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
sat_list = unpack(satellites)
sat_norad_ids = [int(sat.name.split()[1]) for sat in sat_list]

with open(DATA_PATH + FRAME_FILE, "r") as file:
    frames = file.readlines()


found_ids = dict()

for frame in frames:
    if frame.startswith("IRA"):
        sat_id, rel_time, freq, ira_data = decompose_ira_frame(frame)
        beam_id, lat, lon, alt, x, y, z = ira_data

        if not (sat_id and rel_time and lat is not False and lon is not False and alt is not False and alt > 600):
            continue
    else:
        continue

    if not found_ids.get(sat_id):
        found_ids[sat_id] = list()

    # try to identify
    time = Time(START_TIME, format="iso", scale="utc") + TimeDelta(rel_time * unit.ms)
    satrecs = SatrecArray([sat.satrec for sat in sat_list])
    jds = np.array([time.jd1])
    jfs = np.array([time.jd2])

    # do SGP4
    errcodes, positions, velocities = satrecs.sgp4(jds, jfs)
    pos_teme = CartesianRepresentation(x=positions[:, 0, 0], y=positions[:, 0, 1], z=positions[:, 0, 2], unit=unit.km)
    vel_teme = CartesianDifferential(velocities[:, 0, 0], velocities[:, 0, 1], velocities[:, 0, 2], unit=unit.km / unit.s)
    r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=time)
    pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))
    geo_pos = pos_itrs.earth_location.geodetic

    # find closest satellites - first by lon then lat
    found_sat_list = np.array([sat_norad_ids, geo_pos.lon, geo_pos.lat,
                               np.sqrt(np.square(lon * unit.deg - geo_pos.lon) + np.square(lat * unit.deg - geo_pos.lat))], dtype="float").T
    found_sat_list = found_sat_list[found_sat_list[:, 3].argsort()]
    found_ids[sat_id].append(int(found_sat_list[0, 0]))


for k in found_ids.keys():
    out = ""
    ids, counts = np.unique(found_ids[k], return_counts=True)

    print(f"{k: 3d}: " + "\t".join(f"Iridium {sat_id: 3d} x {c}" for sat_id, c in zip(ids, counts)))

for k in found_ids.keys():
    out = ""
    ids, counts = np.unique(found_ids[k], return_counts=True)
    for sat_id, c in zip(ids, counts):
        print(f"{sat_id}: {k: 3d} ({c})")
