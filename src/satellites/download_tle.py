"""
    Download TLEs from CelesTrak and parse them
"""

import requests
import re
from astropy.time import Time, TimeDelta
from astropy import units as unit
from typing import Dict, List

from src.satellites.satellite import Satellite

DOWNLOAD_PATH = "download\\"
DOWNLOAD_EXTENSION = ".txt"
DOWNLOAD_TIME_FILE = "last_download_time"
DOWNLOAD_EXPIRATION = 8  # hours

SATELLITE_CONSTELLATIONS = ("Iridium", "Orbcomm", "Globalstar")

CONSTELLATION_URL = {"Iridium": "iridium-next",
                     "Globalstar": "globalstar",
                     "Orbcomm": "orbcomm"
                     }


def download_tles(constellations=SATELLITE_CONSTELLATIONS, offline_dir=None) -> Dict[str, Dict[str, Satellite]]:
    tles = dict()

    for constellation in constellations:
        tles[constellation] = _download_and_parse_tles(constellation, offline_dir)

    return tles


def unpack(tles: dict) -> List[Satellite]:
    satellite_list = list()
    for const in tles.values():
        satellite_list.extend(const.values())
    return satellite_list


def _download_and_parse_tles(constellation, offline_dir=None):
    print(f"Downloading TLEs for {constellation}...")
    retval = dict()

    tle_lines = _get_tles(constellation, offline_dir)

    for i in range(0, len(tle_lines), 3):
        pattern = re.compile(r"([A-Za-z]+)[- ]([A-Za-z0-9]+)( \[.\])*", re.IGNORECASE)
        match = pattern.match(tle_lines[i].strip())
        if match is None:
            raise ValueError(f"Error in parsing TLEs. Loaded TLEs follow:\n{tle_lines}")

        _, identifier, status = match.groups()

        if status is None:
            status = True
        elif status.lstrip().replace("[", "").replace("]", "") in ["+"]:  # satellite is active
            status = True
        else:
            status = False

        if status:
            retval[identifier] = Satellite(constellation, identifier, tle_lines[i + 1], tle_lines[i + 2])

    print(f"Downloaded")
    return retval


def _get_tles(constellation, offline_dir=None):
    found, expired, timestamp, tles_offline = _read_saved_tles(constellation,
                                                               offline_dir if offline_dir else DOWNLOAD_PATH)

    if found and (not expired or offline_dir):  # up-to-date backup exists
        print("Using previously downloaded TLEs")
        return tles_offline
    else:
        tles = _download_tles_from_server(constellation)

        if tles:  # download ok
            _save_tles_to_file(constellation, tles)
            return tles
        elif found:  # download failed but backup exists
            print(f"There was error in downloading from server, falling back to downloaded version from {timestamp}")
            return tles_offline
        else:  # download failed and no backup exists
            raise RuntimeError("There was error in downloading from server and no previous copy was found")


def _download_tles_from_server(constellation):
    url = f"https://celestrak.org/NORAD/elements/gp.php?GROUP={CONSTELLATION_URL[constellation]}&FORMAT=tle"

    response = requests.get(url)

    if response.status_code >= 300:
        print(f"Error {response.status_code}: {response.reason} while fetching TLEs from {url}")
        return False
    else:
        return response.text.splitlines()


def _save_tles_to_file(constellation, tles):
    with open(DOWNLOAD_PATH + constellation.replace(" ", "_") + DOWNLOAD_EXTENSION, "w+") as file:
        file.write(str(Time.now()) + "\n")
        file.writelines("\n".join(tles))


def _read_saved_tles(constellation, path=DOWNLOAD_PATH):
    try:
        with open(path + constellation.replace(" ", "_") + DOWNLOAD_EXTENSION, "r") as file:
            content = file.readlines()
            timestamp = Time(content[0].strip())
            content = content[1:]
    except FileNotFoundError:
        return False, None, None, None

    # check expiration
    if Time.now() - timestamp > TimeDelta(DOWNLOAD_EXPIRATION * 3600 * unit.s):
        expired = True
    else:
        expired = False

    return True, expired, timestamp, content
