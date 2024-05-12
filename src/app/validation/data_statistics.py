from astropy.coordinates import ITRS
from astropy.time import Time, TimeDelta
from astropy import units as unit
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

from src.utils.data import get_fig_filename
from src.utils.printing import tableify
from src.utils.run_for_all_data import run_for_all_data, load_results
from src.config.locations import LOCATIONS
from src.navigation.calculations import latlon_distance
from src.navigation.curve_fit_method import solve
from src.navigation.data_processing import find_curves
from src.config.parameters import CurveFitMethodParameters
from src.navigation.data_processing import NavDataArrayIndices as IDX
from src.config.setup import *

# first: 5, alt
# second: 5, no alt
# third: 5, alt, no eststate
RUN_NEW_DATA = False
TIME_STEP = 5  # min
DATA_IDX = 0
EXCLUDED_DATA = []

home_lon, home_lat, home_alt = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]
LON_HOME, LAT_HOME, ALT_HOME = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]
plt.style.use("Plots/ctuthesis.mplstyle")


class DataStatistics:
    def __init__(self):
        self.name = None
        self.start_time = None
        self.duration = None
        self.length = None
        self.min_fps = None
        self.mean_fps = None
        self.max_fps = None
        self.received_sats = None
        self.svt_min = None
        self.svt_mean = None
        self.svt_max = None
        self.tlg_min = None
        self.tlg_mean = None
        self.tlg_max = None
        self.min_lon = None
        self.max_lat = None
        self.max_lon = None
        self.min_lat = None
        self.avg_lat = None
        self.avg_lon = None
        self.extent_n = None
        self.extent_s = None
        self.extent_e = None
        self.extent_w = None
        self.cnt_ira = None
        self.cnt_ibc = None
        self.cnt_all = None
        self.mean_noise = None
        self.mean_conf = None

    def print_data(self):
        for k, v in self.__dict__.items():
            print(k, ": ", v)


# section: ------------------------------------------------------------------------- SOLVING
def min_mean_max(arr):
    return np.min(arr), np.mean(arr), np.max(arr)


def cb_solve(nav_data, satellites, default_parameters: CurveFitMethodParameters, exp_name):
    data_path = WORKING_DIR + f"Data\\{exp_name}\\"
    results = DataStatistics()
    results.name = exp_name

    # time ang length
    with open(data_path + "start_time.txt", "r") as file:
        start_time = file.readlines()[0].strip()
    results.start_time = Time(start_time)
    results.duration = nav_data[-1, IDX.t] - nav_data[0, IDX.t]  # s
    results.length = nav_data.shape[0]

    # sats and fps
    sats, frames_per_sat = np.unique(nav_data[:, IDX.sat_id], return_counts=True)
    results.min_fps, results.mean_fps, results.max_fps = min_mean_max(frames_per_sat)
    results.received_sats = ", ".join(str(int(s)) for s in sats)
    sat_vis_times = [nav_data[nav_data[:, IDX.sat_id] == s, IDX.t].max() - nav_data[nav_data[:, IDX.sat_id] == s, IDX.t].min() for s in sats]
    results.svt_min, results.svt_mean, results.svt_max = min_mean_max(sat_vis_times)

    # tle ages
    tle_ages = [satellites[str(int(s))].tle_age(Time(start_time)) for s in sats]
    results.tlg_min, results.tlg_mean, results.tlg_max = min_mean_max(tle_ages)

    # satellite geographic extent
    curves = find_curves(nav_data, min_curve_length=10, max_time_gap=60)
    min_lon, max_lat, max_lon, min_lat = False, False, False, False
    avg_lon, avg_lat = 0, 0
    for curve_array in curves:
        print("#", end="")
        sat_track = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km,
                         z=curve_array[:, IDX.z] * unit.km,
                         v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
                         v_z=curve_array[:, IDX.vz] * unit.km / unit.s).earth_location.geodetic
        min_lon = min(sat_track.lon.min().value, min_lon) if min_lon is not False else sat_track.lon.min().value
        max_lat = max(sat_track.lat.max().value, max_lat) if max_lat is not False else sat_track.lat.max().value
        max_lon = max(sat_track.lon.max().value, max_lon) if max_lon is not False else sat_track.lon.max().value
        min_lat = min(sat_track.lat.min().value, min_lat) if min_lat is not False else sat_track.lat.min().value
        avg_lon += sat_track.lon.mean().value
        avg_lat += sat_track.lat.mean().value
    avg_lat /= len(curves)
    avg_lon /= len(curves)

    results.avg_lat, results.avg_lon = avg_lat, avg_lon
    results.min_lon = min_lon
    results.max_lat = max_lat
    results.max_lon = max_lon
    results.min_lat = min_lat

    results.extent_n = latlon_distance(home_lat, max_lat, home_lon, home_lon)
    results.extent_s = latlon_distance(home_lat, min_lat, home_lon, home_lon)
    results.extent_e = latlon_distance(home_lat, home_lat, home_lon, max_lon)
    results.extent_w = latlon_distance(home_lat, home_lat, home_lon, min_lon)

    # frame statistics
    with open(data_path + "decoded.txt", "r") as file:
        frames = file.readlines()
    cnt_ira, cnt_ibc = 0, 0
    fr_noise, fr_conf = list(), list()
    for f in frames:
        if f.startswith("IRA"):
            cnt_ira += 1
            conf = int(f.split()[4].replace("%", ""))
            noise = [float(n) for n in f.split()[5].split("|")]
            fr_conf.append(conf)
            fr_noise.append(noise)
        elif f.startswith("IBC"):
            cnt_ibc += 1
    results.cnt_ira, results.cnt_ibc = cnt_ira, cnt_ira
    results.cnt_all = len(frames)
    results.mean_noise = np.array(fr_noise).mean(axis=0)
    results.mean_conf = np.mean(fr_conf)

    results.print_data()
    return results


# section: ------------------------------------------------------------------------- ANALYSIS
results = load_results("data_statistics", index=DATA_IDX)
if RUN_NEW_DATA or results is None:
    run_for_all_data(cb_solve, "data_statistics")
    results = load_results("data_statistics")

# section: data
for excl in EXCLUDED_DATA:
    print(f"Excluding {excl}")
    results.pop(excl)

min_results, mean_results, max_results = DataStatistics(), DataStatistics(), DataStatistics()
exps = results.keys()
for res_name in mean_results.__dict__.keys():
    ress = list()
    for exp in exps:
        ress.append(results[exp].__dict__[res_name])
    try:
        min_res = np.min(ress)
        mean_res = np.mean(ress)
        max_res = np.max(ress)

        min_results.__dict__[res_name] = min_res
        mean_results.__dict__[res_name] = mean_res
        max_results.__dict__[res_name] = max_res

    except Exception:
        print(f"Cannot parse {res_name}")
        pass

# table of results
table_data = [r.__dict__.values() for r in results.values()]
tableify(table_data,
         col_head=list(mean_results.__dict__.keys()),
         row_head=list(results.keys()),
         col_dec=[0]*len(list(mean_results.__dict__.keys())))

# table of mean result
table_data = list()
for res_name in mean_results.__dict__.keys():
    table_data.append([min_results.__dict__[res_name], mean_results.__dict__[res_name], max_results.__dict__[res_name], ""])

tableify(table_data,
         row_head=list(mean_results.__dict__.keys()),
         col_head=["", "Min", "Mean", "Max", "Unit"],
         col_dec=[0]*len(table_data))



