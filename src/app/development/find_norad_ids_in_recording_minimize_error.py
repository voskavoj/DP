from astropy.time import Time, TimeDelta
import astropy.units as unit
from sgp4.api import SatrecArray

from src.config.setup import *
from src.radio.iridium_frame_operations import decompose_ira_frame
from src.satellites.download_tle import download_tles, unpack
from src.utils.data import dump_data
from src.satellites.predictions import find_closest_tle_id_to_ira_id_by_position as compare_ira_vs_tle_positions


satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
sat_list = unpack(satellites)
sat_norad_ids = [int(sat.name.split()[1]) for sat in sat_list]

with open(DATA_PATH + FRAME_FILE, "r") as file:
    frames = file.readlines()


found_ids = dict()
ress = list()
base_offset = 0
direction = 0
avg_offset = 0

for frame in frames:
    if frame.startswith("IRA"):
        sat_id, rel_time, freq, ira_data = decompose_ira_frame(frame)
        beam_id, lat, lon, alt, x, y, z = ira_data

        if not (sat_id and rel_time and lat is not False and lon is not False and alt is not False and 700 < alt < 900):
            continue
    else:
        continue

    if not found_ids.get(sat_id):
        found_ids[sat_id] = list()

    time = Time(START_TIME, format="iso", scale="utc") + TimeDelta(rel_time * unit.ms) + TimeDelta(base_offset * unit.s)
    satrecs = SatrecArray([sat.satrec for sat in sat_list])

    # identify
    tle_id, dist, *_ = compare_ira_vs_tle_positions(time, satrecs, sat_list, lat, lon, alt)

    # pick one satellite
    satrecs = SatrecArray([satellites["Iridium"][str(tle_id)].satrec])

    for step in [5, 1, 0.1]:
        time_p = Time(START_TIME, format="iso", scale="utc") + TimeDelta(rel_time * unit.ms) + TimeDelta(base_offset * unit.s) - TimeDelta(step * unit.s)
        _, dist_p, *_ = compare_ira_vs_tle_positions(time_p, satrecs, sat_list, lat, lon, alt)
        time_f = Time(START_TIME, format="iso", scale="utc") + TimeDelta(rel_time * unit.ms) + TimeDelta(base_offset * unit.s) + TimeDelta(step * unit.s)
        _, dist_f, *_ = compare_ira_vs_tle_positions(time_f, satrecs, sat_list, lat, lon, alt)

        if dist <= dist_p and dist <= dist_f:
            direction = 0
        elif dist_p < dist_f:
            direction = -step
        else:
            direction = +step

        while direction != 0:
            time += TimeDelta(direction * unit.s)
            tle_id, dist_1, *_ = compare_ira_vs_tle_positions(time, satrecs, sat_list, lat, lon, alt)

            if dist <= dist_1:
                direction = 0
            else:
                dist = dist_1
                base_offset += direction

    ress.append([sat_id, base_offset, dist])
    avg_offset += base_offset
    print(f"Iri ID {sat_id:03d}, TLE ID {tle_id:03d}, distance {dist / 1000:.0f} km, offset {base_offset:.1f} s (avg {avg_offset/len(ress):.1f} s)")

dump_data("iridium_id_minimization_results", ress)
