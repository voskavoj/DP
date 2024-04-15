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

# ---------------------------- init
satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
sat_list = unpack(satellites)
sat_norad_ids = [int(sat.name.split()[1]) for sat in sat_list]

with open(DATA_PATH + FRAME_FILE, "r") as file:
    frames = file.readlines()


found_ids = dict()
ress = list()
# frames = []
for i, frame in enumerate(frames):
    if i % 500 == 0:
        print(f"\r{100 * i / len(frames):.1f}%", end="")
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

    iri_id = int(sat_id)
    closest_sat_id = int(found_sat_list[0, 0])
    frame_lat, frame_lon, frame_alt = lat, lon, alt
    distances = [latlon_distance(frame_lat, lat.value, frame_lon, lon.value, frame_alt, alt.value) for lat, lon, alt in
                 zip(geo_pos.lat, geo_pos.lon, geo_pos.height)]
    min_idx = np.argmin(np.array(distances))
    closest_sat_id_2 = sat_list[min_idx].number

    dist = distances[min_idx]
    iri_lat, iri_lon, iri_alt = geo_pos.lat[min_idx].value, geo_pos.lon[min_idx].value, geo_pos.height[min_idx].value

    # print(f"Orig ID: {closest_sat_id}, new ID: {closest_sat_id_2}")
    # print(f"ID {iri_id:03d} is Iridium {closest_sat_id_2:03d}, distance {distances[min_idx] / 1000:.0f} km")
    # print(f"Frame coordinates: {frame_lat:.2f} deg, {frame_lon:.2f} deg, {frame_alt} km")
    # print(f"Pred. coordinates: {geo_pos.lat[min_idx]:.2f}, {geo_pos.lon[min_idx]:.2f}, {geo_pos.height[min_idx]:.0f}")
    # print()
        # print(f"{iri_id:03d}\t{closest_sat_id_2:03d}\t{dist/1000:04.0f}\t"
        #       f"{frame_lat-iri_lat:.2f}\t{frame_lon-iri_lon:.2f}\t{frame_alt-iri_alt:.0f}")
        # print(f"{iri_id:03d}\t{closest_sat_id_2:03d}\t{dist/1000:04.0f}\t"
        #       f"{deg_lat_to_m(frame_lat-iri_lat)/1000:.0f}\t{deg_lon_to_m(frame_lon-iri_lon, frame_lon+iri_lon/2)/1000:.0f}\t{frame_alt-iri_alt:.0f}")
    ress.append([iri_id, closest_sat_id_2, dist, frame_lat, frame_lon, frame_alt, iri_lat, iri_lon, iri_alt])

dump_data(f"iridium_id_comparison_results_{EXP_NAME}", ress)
ress = load_data(f"iridium_id_comparison_results_{EXP_NAME}_0")

for k in found_ids.keys():
    out = ""
    ids, counts = np.unique(found_ids[k], return_counts=True)

    print(f"{k: 3d}: " + "\t".join(f"Iridium {sat_id: 3d} x {c}" for sat_id, c in zip(ids, counts)))

for k in sorted(found_ids.keys()):
    out = ""
    ids, counts = np.unique(found_ids[k], return_counts=True)
    for sat_id, c in zip(ids, counts):
        print(f"{sat_id}: {k: 3d} ({c})")


import numpy as np
from collections import defaultdict

# Assuming your ress_arr is a numpy array
# Sample data:
ress_arr = np.array(ress)

# Dictionary to store results
result = defaultdict(lambda: defaultdict(list))

# Iterate over rows of ress_arr
for row in ress_arr:
    iri_id, closest_sat_id_2, dist = row[0], row[1], row[2]
    result[closest_sat_id_2][iri_id].append(dist)

# Sort the result by occurrences
sorted_result = {}
for sat_id, iri_data in result.items():
    sorted_iri_data = sorted(iri_data.items(), key=lambda x: len(x[1]), reverse=True)
    sorted_result[sat_id] = sorted_iri_data

# Print the sorted result
for sat_id, iri_data in sorted(sorted_result.items()):
    print(f"Iridium {sat_id:.0f}:")
    for iri_id, dist_list in iri_data:
        occurrences = len(dist_list)
        avg_dist = np.mean(dist_list)
        if occurrences > 10 and avg_dist < 100e3:
            dict_style = f" | {sat_id:.0f}: {iri_id:.0f},  # {occurrences}"
        else:
            dict_style = ""
        print(f"  IRI ID: {iri_id: 4.0f}, Occurrences: {occurrences: 4d}, Average Distance: {avg_dist/1000: 7.1f} km" + dict_style)

