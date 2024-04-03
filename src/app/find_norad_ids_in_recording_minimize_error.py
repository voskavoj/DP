import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, ITRS, TEME, EarthLocation
from astropy.time import Time, TimeDelta
import astropy.units as unit
from sgp4.api import SatrecArray

from src.config.setup import *
from src.navigation.calculations import latlon_distance, deg_lat_to_m, deg_lon_to_m
from src.radio.iridium_frame_operations import decompose_ira_frame
from src.satellites.download_tle import download_tles, unpack
from src.utils.data import dump_data, load_data


def compare_ira_vs_tle_positions(time, satrecs, sat_list, iri_id, iri_lat, iri_lon, iri_alt):

    jds = np.array([time.jd1])
    jfs = np.array([time.jd2])
    errcodes, positions, velocities = satrecs.sgp4(jds, jfs)
    pos_teme = CartesianRepresentation(x=positions[:, 0, 0], y=positions[:, 0, 1], z=positions[:, 0, 2], unit=unit.km)
    vel_teme = CartesianDifferential(velocities[:, 0, 0], velocities[:, 0, 1], velocities[:, 0, 2],
                                     unit=unit.km / unit.s)
    r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=time)
    pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))
    geo_pos = pos_itrs.earth_location.geodetic

    distances = [latlon_distance(iri_lat, lat.value, iri_lon, lon.value, iri_alt, alt.value) for lat, lon, alt in
                 zip(geo_pos.lat, geo_pos.lon, geo_pos.height)]
    min_idx = np.argmin(np.array(distances))
    closest_sat_id_2 = sat_list[min_idx].number

    dist = distances[min_idx]
    tle_lat, tle_lon, tle_alt = geo_pos.lat[min_idx].value, geo_pos.lon[min_idx].value, geo_pos.height[min_idx].value

    return iri_id, closest_sat_id_2, dist, iri_lat, iri_lon, iri_alt, tle_lat, tle_lon, tle_alt


satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
sat_list = unpack(satellites)
sat_norad_ids = [int(sat.name.split()[1]) for sat in sat_list]

with open(DATA_PATH + FRAME_FILE, "r") as file:
    frames = file.readlines()


found_ids = dict()
ress = list()
base_offset = 0
direction = 0

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

    time = Time(START_TIME, format="iso", scale="utc") + TimeDelta(rel_time * unit.ms) + TimeDelta(base_offset * unit.s)
    satrecs = SatrecArray([sat.satrec for sat in sat_list])

    # identify
    _, tle_id, dist, *_ = compare_ira_vs_tle_positions(time, satrecs, sat_list, int(sat_id), lat, lon, alt)

    # pick one satellite
    satrecs = SatrecArray([satellites["Iridium"][str(tle_id)].satrec])

    for step in [5, 1, 0.1]:
        time_p = Time(START_TIME, format="iso", scale="utc") + TimeDelta(rel_time * unit.ms) + TimeDelta(base_offset * unit.s) - TimeDelta(step * unit.s)
        _, _, dist_p, *_ = compare_ira_vs_tle_positions(time_p, satrecs, sat_list, int(sat_id), lat, lon, alt)
        time_f = Time(START_TIME, format="iso", scale="utc") + TimeDelta(rel_time * unit.ms) + TimeDelta(base_offset * unit.s) + TimeDelta(step * unit.s)
        _, _, dist_f, *_ = compare_ira_vs_tle_positions(time_f, satrecs, sat_list, int(sat_id), lat, lon, alt)

        if dist <= dist_p and dist <= dist_f:
            direction = 0
        elif dist_p < dist_f:
            direction = -step
        else:
            direction = +step

        while direction != 0:
            time += TimeDelta(direction * unit.s)
            _, tle_id, dist_1, *_ = compare_ira_vs_tle_positions(time, satrecs, sat_list, int(sat_id), lat, lon, alt)

            if dist <= dist_1:
                direction = 0
            else:
                dist = dist_1
                base_offset += direction

    print(f"Iri ID {sat_id:03d}, TLE ID {tle_id:03d}, distance {dist / 1000:.0f} km, offset {base_offset:.1f} s")
    ress.append([sat_id, base_offset, dist])

dump_data("iridium_id_minimization_results", ress)







