import glob
import os

from astropy.time import Time

from src.satellites.satellite import Satellite, IRIDIUM_NORAD_IDS
from src.config.setup import *
from src.satellites.download_tle import CONSTELLATION_URL


def load_tles_from_folder(folder_path: str, start_time: str, allow_newer_tles=False):
    start_time = Time(start_time, format="iso", scale="utc")
    tles = dict()
    for filename in glob.glob(os.path.join(folder_path, '*.txt')):
        with open(os.path.join(os.getcwd(), filename), 'r') as tle_file:
            tle_lines = tle_file.readlines()
            for i, line in enumerate(tle_lines):
                if line.startswith("1 "):
                    tmp_sat = Satellite("", "0", tle_lines[i], tle_lines[i + 1])
                    norad_id = tmp_sat.satrec.satnum
                    if tle_id := IRIDIUM_NORAD_IDS.get(norad_id, False):
                        constellation = "Iridium"
                        # print(tle_lines[i + 1])
                        sat = Satellite(constellation, str(tle_id), tle_lines[i], tle_lines[i + 1])
                    else:
                        continue

                    if constellation not in tles:
                        tles[constellation] = dict()
                    # a TLE set already exists and is less up to date
                    if str(tle_id) in tles[constellation] and abs(sat.tle_age(start_time)) < abs(tles[constellation][str(tle_id)].tle_age(start_time)):
                        print(f"Found TLE set for {constellation} {tle_id} with more up to date set: old age {tles[constellation][str(tle_id)].tle_age(start_time) / 3600:.1f} h, new age {sat.tle_age(start_time) / 3600:.1f} h", end=", ")
                        if allow_newer_tles:
                            print("overwriting")
                            tles[constellation][str(tle_id)] = sat
                        else:
                            if sat.tle_age(start_time) >= 0:
                                print("overwriting")
                                tles[constellation][str(tle_id)] = sat
                            else:
                                print(f"but this is newer than start time - keeping the old")
                    # no TLE set exists
                    elif str(tle_id) not in tles[constellation]:
                        tles[constellation][str(tle_id)] = sat
                    # a TLE set already exists and is more up to date
                    else:
                        pass

    return tles


satellites = load_tles_from_folder(DATA_PATH + "tles\\", START_TIME)

for const in satellites.keys():
    with open(DATA_PATH + "tles\\" + f"{CONSTELLATION_URL[const]}.txt", "w") as file:
        for sat in satellites[const].values():
            file.writelines([sat.name.upper() + "\n", sat.tle_line_1, sat.tle_line_2])
