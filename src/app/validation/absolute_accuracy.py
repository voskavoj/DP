from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

from src.utils.data import get_fig_filename
from src.utils.printing import tableify
from src.utils.run_for_all_data import run_for_all_data, load_results
from src.config.locations import LOCATIONS
from src.navigation.calculations import latlon_distance, m_to_deg_lon, m_to_deg_lat
from src.navigation.curve_fit_method import solve


RUN_NEW_DATA = False

home_lon, home_lat, home_alt = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]
LON_HOME, LAT_HOME, ALT_HOME = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]
plt.style.use("Plots/ctuthesis.mplstyle")


# section: ------------------------------------------------------------------------- SOLVING
def cb_solve(nav_data, satellites, default_parameters, *_):
    res = solve(nav_data, satellites, default_parameters)
    return res


# section: ------------------------------------------------------------------------- ANALYSIS
results = load_results("absolute_accuracy")
if RUN_NEW_DATA or results is None:
    run_for_all_data(cb_solve, "absolute_accuracy")
    results = load_results("absolute_accuracy")

# section: data
for exp_name, res in results.items():
    lat, lon, alt, off, dft = res
    hor_dist = latlon_distance(LAT_HOME, lat, LON_HOME, lon)
    abs_dist = latlon_distance(LAT_HOME, lat, LON_HOME, lon, ALT_HOME, alt)
    #                    0    1    2    3    4    5         6
    results[exp_name] = [lat, lon, alt, off, dft, hor_dist, abs_dist]

    print(f"{exp_name}: {results[exp_name]}")

res_arr = np.array(list(results.values()))
tableify(list(results.values()),
         col_dec=[3, 3, 0, 0, 3, 0, 0],
         col_head=["ID", "Lat", "Lon", "Alt (m)", "Offset (Hz)", "Drift (Hz/s)", "Hor. error (m)", "Abs. error (m)"],
         row_head=list(results.keys()))


# section: error and spread
def get_pos_circle(res_arr, percentile):
    mask = res_arr[:, 6] <= np.percentile(res_arr[:, 6], percentile)
    lat_arr, lon_arr, dist_arr = res_arr[mask, 0], res_arr[mask, 1], res_arr[mask, 6]

    mean_lat = lat_arr.mean()
    mean_lon = lon_arr.mean()

    rel_dist_arr = np.array([latlon_distance(mean_lat, lat, mean_lon, lon) for lat, lon in zip(lat_arr, lon_arr)])
    spread = latlon_distance(mean_lat, lat_arr[rel_dist_arr.argmax()], mean_lon, lon_arr[rel_dist_arr.argmax()])
    err = latlon_distance(mean_lat, LAT_HOME, mean_lon, LON_HOME)

    return mean_lat, mean_lon, err, spread


mean_lat_all, mean_lon_all, err_all, spread_all = get_pos_circle(res_arr, 100)
mean_lat_cep, mean_lon_cep, err_cep, spread_cep = get_pos_circle(res_arr, 50)
mean_lat_95,  mean_lon_95,  err_95,  spread_95 = get_pos_circle(res_arr, 95)

tableify([get_pos_circle(res_arr, 100), get_pos_circle(res_arr, 50), get_pos_circle(res_arr, 95)],
         col_dec=[3, 3, 0, 0], col_head=["", "Mean lat", "Mean lon", "Error (m)", "Spread (m)"], row_head=["All", "CEP", "95%"])


# section: probability
def get_cumulative_distribution_fcn(dist_arr, dist_max, percentile=100):
    mask = dist_arr <= np.percentile(dist_arr, percentile)
    dist_arr = dist_arr[mask]
    dist_bins = (list(dist_arr))
    dist_bins.extend([0, dist_max])
    cumm_func_list = list()
    for db in sorted(dist_bins):
        mask = dist_arr <= db
        good_results = dist_arr[mask]
        prob = len(good_results) / len(dist_arr)
        cumm_func_list.append([db, prob])
    return np.array(cumm_func_list)


def get_cumulative_distribution_fcn_spread(res_arr, dist_max, percentile=100):
    mean_lat, mean_lon, err, spread = get_pos_circle(res_arr, percentile)
    mask = res_arr[:, 6] <= np.percentile(res_arr[:, 6], percentile)
    dist_arr = np.array([latlon_distance(mean_lat, lat, mean_lon, lon) for
                         lat, lon in zip(res_arr[mask, 0], res_arr[mask, 1])])
    return get_cumulative_distribution_fcn(dist_arr.T, dist_max)


# cumulative function for all
cumm_func_all = get_cumulative_distribution_fcn(res_arr[:, 6], 3500)
cumm_func_all_spread = get_cumulative_distribution_fcn_spread(res_arr, 3500, 100)


# cumulative function for CEP
cumm_func_cep = get_cumulative_distribution_fcn(res_arr[:, 6], 3500, 50)
cumm_func_cep_spread = get_cumulative_distribution_fcn_spread(res_arr, 3500, 50)

# cumulative function for CEP
cumm_func_95 = get_cumulative_distribution_fcn(res_arr[:, 6], 3500, 95)
cumm_func_95_spread = get_cumulative_distribution_fcn_spread(res_arr, 3500, 95)


# section: plotting
def make_circle(center_x, center_y, radius):
    point_count = int(radius / 1000 * 400)
    angles = np.linspace(0, 2 * np.pi, point_count)
    x = center_x + radius * np.cos(angles)
    y = center_y + radius * np.sin(angles)
    return x, y


def make_latlon_circle(center_lat, center_lon, radius):
    point_count = int(radius / 1000 * 400)
    angles = np.linspace(0, 2 * np.pi, point_count)
    x = center_lon + m_to_deg_lon(radius, center_lat) * np.cos(angles)
    y = center_lat + m_to_deg_lat(radius, center_lat) * np.sin(angles)
    return x, y


plt.figure()
plt.figtext(.12, .92, f"Grid size is 1 km")
ll_lon = LON_HOME - m_to_deg_lon(4e3, home_lat)
ll_lat = LAT_HOME - m_to_deg_lat(4e3, home_lat)
ur_lon = LON_HOME + m_to_deg_lon(4e3, home_lat)
ur_lat = LAT_HOME + m_to_deg_lat(4e3, home_lat)

m = Basemap(llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat,
            rsphere=(6378137.00, 6356752.3142), resolution='h', projection='merc', )
# positions
m.plot(home_lon, home_lat,           "x", color="blue",   latlon=True, label="Actual position")
m.plot(res_arr[:, 1], res_arr[:, 0], "o", color="orange", latlon=True, label="Estimated positions")
# spread
m.plot(mean_lon_all, mean_lat_all, marker=".", color="red",   latlon=True, label="Mean all")
m.plot(mean_lon_cep, mean_lat_cep, marker=".", color="green", latlon=True, label="Mean CEP")
m.plot(mean_lon_95,  mean_lat_95,  marker=".", color="blue",  latlon=True, label="Mean 95%")
m.plot(*make_latlon_circle(mean_lat_all, mean_lon_all, spread_all), ".", ms=1, color="red",   label="_nolegend_", alpha=0.5, latlon=True)
m.plot(*make_latlon_circle(mean_lat_cep, mean_lon_cep, spread_cep), ".", ms=1, color="green", label="_nolegend_", alpha=0.5, latlon=True)
m.plot(*make_latlon_circle(mean_lat_95,  mean_lon_95,  spread_95),  ".", ms=1, color="blue",  label="_nolegend_", alpha=0.5, latlon=True)
# distance around home
m.plot(*make_latlon_circle(home_lat, home_lon, 1000), ".", ms=1, color="black", label="_nolegend_", alpha=0.25, latlon=True)
m.plot(*make_latlon_circle(home_lat, home_lon, 2000), ".", ms=1, color="black", label="_nolegend_", alpha=0.25, latlon=True)
m.plot(*make_latlon_circle(home_lat, home_lon, 3000), ".", ms=1, color="black", label="_nolegend_", alpha=0.25, latlon=True)
m.drawcoastlines()
m.drawcountries()
m.drawrivers(color="blue")
m.drawparallels(np.arange(ll_lat, ur_lat, m_to_deg_lat(1e3, home_lat)), labels=[0, 1, 0, 0])
m.drawmeridians(np.arange(ll_lon, ll_lat, m_to_deg_lon(1e3, home_lat)), labels=[0, 0, 0, 1], rotation="vertical")
plt.legend()
plt.tight_layout()
plt.savefig(get_fig_filename("validation\\" + f"absolute_accuracy", idx=False), dpi=600)


# plt.figure()
# plt.figtext(.15, .05, f"")
# plt.plot(deg_lon_to_m(res_arr[:, 1] - home_lon, home_lat), deg_lat_to_m(res_arr[:, 0] - home_lat, home_lat), marker=".", color="orange", label="Estimated positions")
# plt.plot(0, 0, "x", color="blue", label="Actual position")
# plt.plot(*make_circle(0, 0, 1000), ".", ms=1, color="black", label="_nolegend_", alpha=0.25)
# plt.plot(*make_circle(0, 0, 2000), ".", ms=1, color="black", label="_nolegend_", alpha=0.25)
# plt.plot(*make_circle(0, 0, 3000), ".", ms=1, color="black", label="_nolegend_", alpha=0.25)
# plt.grid()
# plt.axis('square')
# plt.legend()

plt.figure()
plt.plot(cumm_func_all[:, 0], cumm_func_all[:, 1], ".-", label="To actual postion")
plt.plot(cumm_func_all_spread[:, 0], cumm_func_all_spread[:, 1], ".-", label="To mean postion")
plt.xlabel("Error (m)")
plt.ylabel("Probability (-)")
plt.legend()
plt.savefig(get_fig_filename("validation\\" + f"absolute_accuracy_all", idx=False), dpi=600)

plt.figure()
plt.plot(cumm_func_cep[:, 0],        cumm_func_cep[:, 1], ".-", label="To actual postion")
plt.plot(cumm_func_cep_spread[:, 0], cumm_func_cep_spread[:, 1], ".-", label="To mean postion")
plt.xlabel("Error (m)")
plt.ylabel("Probability (-)")
plt.legend()
plt.savefig(get_fig_filename("validation\\" + f"absolute_accuracy_cep", idx=False), dpi=600)

plt.figure()
plt.plot(cumm_func_95[:, 0],        cumm_func_95[:, 1], ".-", label="To actual postion")
plt.plot(cumm_func_95_spread[:, 0], cumm_func_95_spread[:, 1], ".-", label="To mean postion")
plt.xlabel("Error (m)")
plt.ylabel("Probability (-)")
plt.legend()
plt.savefig(get_fig_filename("validation\\" + f"absolute_accuracy_95", idx=False), dpi=600)


plt.show()
