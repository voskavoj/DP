from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

from src.navigation.calculations import latlon_distance
from src.navigation.data_processing import nav_data_to_array
from src.utils.data import load_data, get_fig_filename
from src.config.locations import LOCATIONS


def plot_results_of_iterative_position_finding(data: str | list, r=None, show=False):
    if isinstance(data, str):
        results = load_data(data)
    else:
        results = data

    final_lat, final_lon = min(results, key=lambda x: x[0])[2], min(results, key=lambda x: x[0])[3]

    home_lon, home_lat = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1]

    print(f"Position error: {latlon_distance(home_lat, final_lat, home_lon, final_lon):.1f} m")
    
    plt.figure()
    m = Basemap(llcrnrlon=0, llcrnrlat=40, urcrnrlon=30, urcrnrlat=65,
                rsphere=(6378137.00, 6356752.3142), resolution='h', projection='merc', )
    
    res_arr = np.array(results)
    m.plot(home_lon, home_lat, latlon=True, marker="x", label="Actual position")
    m.plot(final_lon, final_lat, latlon=True, marker="o", label="Estimated position")
    m.plot(res_arr[:, 3], res_arr[:, 2], latlon=True, label="Algorithm path")
    if r is not None:
        sat_track = r.earth_location.geodetic
        m.plot(sat_track.lon, sat_track.lat, latlon=True, label="Satellite track")
    m.drawcoastlines()
    m.fillcontinents()
    m.drawcountries()
    m.drawrivers(color="blue")
    m.drawparallels(np.arange(40, 65, 5), labels=[0, 1, 0, 0])
    m.drawmeridians(np.arange(0, 30, 5), labels=[0, 0, 0, 1])
    plt.legend()

    if show:
        plt.show()


def plot_analyzed_curve(curve, dopp_start, dopp_end, curve_duration, curve_density, largest_gap, variance):
    curve_array = nav_data_to_array(curve)
    plt.figure()
    plt.title(f"T:{(dopp_start > 0 > dopp_end or dopp_start < 0 < dopp_end)}, "
              f"L:{curve_duration:.1f}, "
              f"D:{curve_density:.3f} \n"
              f"G{largest_gap:.2f}, "
              f"V{variance:.3f}\n")
    plt.plot(curve_array[:, 0], curve_array[:, 1] - curve_array[:, 2], ".")
    plt.savefig(get_fig_filename("analyzed_curve"))
    plt.close()
