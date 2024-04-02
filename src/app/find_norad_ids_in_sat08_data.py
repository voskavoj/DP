import re
import glob
import numpy as np

from astropy.coordinates import CartesianRepresentation, CartesianDifferential, ITRS, TEME
from astropy.time import Time, TimeDelta
import astropy.units as unit
from sgp4.api import SatrecArray

from src.navigation.calculations import latlon_distance
from src.satellites.satellite import Satellite

DATA_FILE_PATH = "Data\\sat08\\1109-0910_20_parsed.txt"
SAVE_FILE_PATH = "Data\\sat08\\iridium_ids.csv"
TLE_PATH = "Data\\sat08\\tle\\"
IRIDIUM_IDS = [2, 3, 4, 5, 6, 7, 8, 9, 13, 16, 17, 18, 22, 23, 24, 25, 26, 28, 29, 30, 33, 36, 38, 39, 40, 42, 43, 44, 46, 48, 49, 50, 51, 57, 65, 67, 68, 69, 71, 72, 73, 74, 77, 78, 79, 81, 82, 85, 87, 88, 89, 90, 92, 93, 94, 96, 99, 103, 104, 107, 109, 110, 111, 112, 114, 115]


# Load TLEs from Special Request
satellites_all = dict()
for tle_file in glob.glob(TLE_PATH + "*.txt"):
    with open(tle_file, "r") as file:
        tle_lines = file.readlines()

        if len(tle_lines) < 3 and not tle_lines[0].upper().startswith("IRIDIUM"):
            # print(f"Skipping {tle_file}")
            continue
        for i in range(0, len(tle_lines), 3):
            pattern = re.compile(r"([A-Za-z]+)[- ]([A-Za-z0-9]+)( \[.\])*", re.IGNORECASE)
            match = pattern.match(tle_lines[i].strip())
            if match is None:
                raise ValueError(f"Error in parsing TLEs. Loaded TLEs follow:\n{tle_lines}")

            _, identifier, status = match.groups()

            if identifier not in satellites_all:
                satellites_all[identifier] = list()
            satellites_all[identifier].append(Satellite("Iridium", identifier, tle_lines[i + 1], tle_lines[i + 2]))
print("TLEs loaded.")

# Create ID table
id_table = np.zeros((200, 120), dtype="int")
last_time_tles_were_chosen = None

# Load data file
with open(DATA_FILE_PATH, "r") as file:
    for line_num, line in enumerate(file):
        if line_num < 0:  # skip already processed lines
            continue

        try:
            # parse line
            # line: timestamp unix | time ms | sat ID | beam ID | lat | lon | alt | stuff...
            time_unix, time_ms, iri_id, _, frame_lat, frame_lon, frame_alt, *_ = line.split()
            time_unix, time_ms, iri_id, frame_lat, frame_lon, frame_alt = int(time_unix), int(time_ms), int(iri_id), float(frame_lat), float(frame_lon), int(frame_alt)
            frame_time = Time(time_unix, format="unix", scale="utc") + TimeDelta(int(time_ms) * unit.ms)

            if not int(iri_id) in IRIDIUM_IDS:
                print("Unknown ID: ", iri_id)
                continue

            if frame_alt < 700:
                continue

            if last_time_tles_were_chosen is None or frame_time - last_time_tles_were_chosen > 5 * unit.hour:
                print("Selecting TLEs")
                last_time_tles_were_chosen = frame_time
                sat_list = list()
                # get TLEs for this frame time
                for sat_id, tle_list in satellites_all.items():
                    # pick from sat_list the closest in time
                    time_diff = np.array([abs((sat.tle_time - frame_time).to_value(unit.s)) for sat in tle_list])
                    sat_list.append(tle_list[np.argmin(time_diff)])

            # do SGP4
            satrecs = SatrecArray([sat.satrec for sat in sat_list])
            jds = np.array([frame_time.jd1])
            jfs = np.array([frame_time.jd2])
            errcodes, positions, velocities = satrecs.sgp4(jds, jfs)
            pos_teme = CartesianRepresentation(x=positions[:, 0, 0], y=positions[:, 0, 1], z=positions[:, 0, 2],
                                               unit=unit.km)
            vel_teme = CartesianDifferential(velocities[:, 0, 0], velocities[:, 0, 1], velocities[:, 0, 2],
                                             unit=unit.km / unit.s)
            r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=frame_time)
            pos_itrs = r_teme.transform_to(ITRS(obstime=frame_time))
            geo_pos = pos_itrs.earth_location.geodetic

            # find the closest match
            distances = [latlon_distance(frame_lat, lat.value, frame_lon, lon.value, frame_alt, alt.value) for lat, lon, alt in zip(geo_pos.lat, geo_pos.lon, geo_pos.height)]
            min_idx = np.argmin(np.array(distances))
            closest_sat_id = sat_list[min_idx].number

            print(f"ID {iri_id:03d} is Iridium {closest_sat_id:03d}, distance {distances[min_idx]/1000:.0f} km, (line {line_num})")
            # print(f"Frame coordinates: {frame_lat:.2f} deg, {frame_lon:.2f} deg, {frame_alt} km")
            # print(f"Pred. coordinates: {geo_pos.lat[min_idx]:.2f}, {geo_pos.lon[min_idx]:.2f}, {geo_pos.height[min_idx]:.0f}")
            # print()

            id_table[closest_sat_id, iri_id] += 1

        except KeyboardInterrupt:
            break

np.savetxt(SAVE_FILE_PATH, id_table, delimiter=";", fmt="%d")
