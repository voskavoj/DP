"""
    Plotting helpers to reduce clutter
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

from src.navigation.calculations import latlon_distance
from src.navigation.data_processing import nav_data_to_array
from src.utils.data import load_data, get_fig_filename
from src.config.locations import LOCATIONS


class ColorPicker:
    def_list = ["b", "g", "r", "c", "m", "y", "k", "aqua", "aquamarine", "blue", "blueviolet",
                "brown", "burlywood", "cadetblue", "chartreuse", "chocolate",
                "coral", "cornflowerblue", "crimson", "cyan", "darkblue", "darkcyan", "darkgoldenrod", "darkgray",
                "darkgreen", "darkgrey", "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", "darkorchid",
                "darkred", "darksalmon", "darkseagreen", "darkslateblue", "darkslategray", "darkslategrey",
                "darkturquoise", "darkviolet", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue",
                "firebrick", "forestgreen", "fuchsia", "gold", "goldenrod", "gray", "green", "greenyellow",
                "grey", "hotpink", "indianred", "indigo", "khaki", "lavender", "lavenderblush", "lawngreen",
                "lightblue", "lightcoral", "lightcyan", "lightgreen", "lightpink", "lightsalmon", "lightseagreen",
                "lightskyblue", "lightslategray", "lightslategrey", "lightsteelblue", "lightyellow", "lime",
                "limegreen", "magenta", "maroon", "mediumaquamarine", "mediumblue", "mediumorchid", "mediumpurple",
                "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred",
                "midnightblue", "mintcream", "mistyrose", "moccasin", "navajowhite", "navy", "oldlace", "olive",
                "olivedrab", "orange", "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise",
                "palevioletred", "peachpuff", "peru", "pink", "plum", "powderblue", "purple", "rebeccapurple", "red",
                "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown", "seagreen", "sienna", "silver",
                "skyblue", "slateblue", "slategray", "slategrey", "springgreen", "steelblue", "tan", "teal",
                "thistle", "tomato", "turquoise", "violet", "yellow", "yellowgreen"]

    def __init__(self):
        self.colors_by_line = dict()
        self.idx = 0
        self.color_list = self.def_list

    def _choose_new_color(self):
        if self.idx >= len(self.color_list):
            print("Warning: recycling colors")
            self.index = 0

        color = self.color_list[self.idx]
        self.idx += 1
        return color

    def reserve_colors(self, *colors):
        for color in colors:
            self.color_list.remove(color)

    def pick_color(self, line_name):
        if line_name in self.colors_by_line:
            color = self.colors_by_line[line_name]
        else:
            color = self._choose_new_color()
            print(f"Color for {line_name}: {color}")
            self.colors_by_line[line_name] = color

        return color


def plot_results_of_iterative_position_finding(data: str | list, r=None, show=False, map_covers_tracks=False):
    if isinstance(data, str):
        results = load_data(data)
    else:
        results = data

    plt.style.use("Plots/ctuthesis.mplstyle")

    final_lat, final_lon, final_alt = results[-1][2], results[-1][3], results[-1][4]
    try:
        final_off, final_dft = results[-1][5], results[-1][6]
    except IndexError:
        final_off, final_dft = None, None
    home_lon, home_lat, home_alt = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]
    pos_error = latlon_distance(home_lat, final_lat, home_lon, final_lon, home_alt, final_alt)
    res_arr = np.array(results)

    min_lat = 51
    min_lon = 12
    max_lat = 54
    max_lon = 18

    if map_covers_tracks and r is not None:
        sat_track = r.earth_location.geodetic
        min_lon = min(sat_track.lon.min().value, min_lon)
        max_lat = max(sat_track.lat.max().value, max_lat)
        max_lon = max(sat_track.lon.max().value, max_lon)
        min_lat = min(sat_track.lat.min().value, min_lat)

    plt.figure()
    plt.title(f"Pos. error: {pos_error:.1f} m")
    # plt.figtext(.15, .05, f"Lat.: {final_lat:.2f}°, Lon.: {final_lon:.2f}°, Alt.: {final_alt:.0f} m" +
    #             f", Offset: {final_off / 1e3:.2f} kHz, Drift: {final_dft * 1e3:.0f} mHz/s" if final_off is not None else "")

    m = Basemap(llcrnrlon=min_lon - 2, llcrnrlat=min_lat - 2, urcrnrlon=max_lon + 2, urcrnrlat=max_lat + 2,
                rsphere=(6378137.00, 6356752.3142), resolution='h', projection='merc', )

    if r is not None:
        sat_track = r.earth_location.geodetic
        m.plot(sat_track.lon, sat_track.lat, ".", color="red", ms=1.5, latlon=True, label="Satellite track")
    m.plot(res_arr[:, 3], res_arr[:, 2], marker=".", color="green", latlon=True, label="Algorithm path")
    m.plot(final_lon, final_lat, "o", color="orange", latlon=True, label="Estimated position")
    m.plot(home_lon, home_lat, "x", color="blue", latlon=True, label="Actual position")
    m.drawcoastlines()
    m.fillcontinents()
    m.drawcountries()
    m.drawrivers(color="blue")
    m.drawparallels(np.arange(40, 65, 2), labels=[0, 1, 0, 0])
    m.drawmeridians(np.arange(0, 30, 2), labels=[0, 0, 0, 1])
    plt.legend(loc="upper right")
    plt.tight_layout()

    plt.savefig(get_fig_filename("fig"))

    if show:
        plt.show()


def plot_analyzed_curve(curve, dopp_start, dopp_end, curve_duration, curve_density, largest_gap, variance, ok=None):
    if ok is True:
        color = "green"
    elif ok is False:
        color = "red"
    else:
        color = "blue"

    curve_array = nav_data_to_array(curve)
    plt.figure()
    plt.title(f"T:{(dopp_start > 0 > dopp_end or dopp_start < 0 < dopp_end)}, "
              f"L:{curve_duration:.1f}, "
              f"D:{curve_density:.3f} \n"
              f"G{largest_gap:.2f}, "
              f"V{variance:.3f}\n")
    plt.plot(curve_array[:, 0], curve_array[:, 1] - curve_array[:, 2], ".", color=color)
    plt.savefig(get_fig_filename("analyzed_curve"))
    plt.close()


def plot_measured_vs_trial_curve(measured_curve, trial_curve, lat, lon, alt, off):
    plt.style.use("Plots/ctuthesis.mplstyle")
    plt.figure()
    plt.xlabel("Unix time [s]")
    plt.ylabel("Doppler shift [Hz]")
    # plt.figtext(.15, .12, f"Lat.: {lat:.2f}°, Lon.: {lon:.2f}°, Alt.: {alt:.0f} m, Offset: {off / 1e3:.0f} kHz")
    plt.scatter(measured_curve[:, 0], measured_curve[:, 1], marker=".", label="Measured curve")
    plt.scatter(trial_curve[:, 0], trial_curve[:, 1], marker=".", label="Trial curve")
    plt.legend()


if __name__ == "__main__":
    for i in [1104]:
        plot_results_of_iterative_position_finding(f"results_{i}")
    plt.show()
