import pickle

import numpy as np
from astropy.coordinates import ITRS, EarthLocation
import astropy.units as unit
from matplotlib import pyplot as plt

from src.config.locations import LOCATIONS
from src.config.setup import *
from src.navigation.data_processing import find_curves, NavDataArrayIndices as IDX
from src.satellites.download_tle import download_tles
from src.utils.plots import plot_measured_vs_trial_curve


def check_trial_curve(lat, lon, alt, off, measured_curve, time_arr, r_sat_arr, v_sat_arr, base_freq, plot=False):
    curve_len = measured_curve.shape[0]

    trial_curve = np.empty((curve_len, 2))
    trial_curve[:, 0] = measured_curve[:, 0]  # times are the same
    r_user_arr = (EarthLocation.from_geodetic(lon, lat, alt)
                  .get_itrs().cartesian.without_differentials())

    vs, rs, ru = v_sat_arr, r_sat_arr, r_user_arr * np.ones(curve_len)
    vs = np.array([vs.d_x.value, vs.d_y.value, vs.d_z.value]) * 1000
    rs = np.array([rs.x.value, rs.y.value, rs.z.value]) * 1000
    ru = np.array([ru.x.value, ru.y.value, ru.z.value]) * 1
    f_b = measured_curve[:, 2]

    # Calculate range rate: Ro_dot = (V_s - V_u) * (r_s - r_u) / ||r_s - r_u||
    rel_vel = np.sum(vs * (rs - ru) / np.linalg.norm(rs - ru, axis=0), axis=0)
    # Calculate doppler frequency
    C = 299792458  # m/s
    f_d = -1 * rel_vel * f_b / C

    # Calculate drift
    drf = 1  # Hz/s
    f_drift = (trial_curve[:, 0] - np.min(trial_curve[:, 0])) * drf

    # Construct trial curve
    trial_curve[:, 1] = f_d + off #+ f_drift

    # Calculate variance
    sum_of_squares = np.sum((measured_curve[:, 1] - trial_curve[:, 1]) ** 2)

    if plot:
        plot_measured_vs_trial_curve(measured_curve, trial_curve, lat, lon, alt, off)

    return trial_curve


def find_offset_and_drift_in_data():
    satellites = download_tles(constellations=CONSTELLATIONS, offline_dir=DATA_PATH)
    lon, lat, alt = LOCATIONS[LOCATION][0], LOCATIONS[LOCATION][1], LOCATIONS[LOCATION][2]


    with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
        nav_data = pickle.load(file)
    detected_curves = find_curves(nav_data)
    curve_array = np.vstack(detected_curves)

    r = ITRS(x=curve_array[:, IDX.x] * unit.km, y=curve_array[:, IDX.y] * unit.km, z=curve_array[:, IDX.z] * unit.km,
             v_x=curve_array[:, IDX.vx] * unit.km / unit.s, v_y=curve_array[:, IDX.vy] * unit.km / unit.s,
             v_z=curve_array[:, IDX.vz] * unit.km / unit.s)
    measured_curve = np.column_stack(
        (curve_array[:, IDX.t], curve_array[:, IDX.f] - curve_array[:, IDX.fb], curve_array[:, IDX.fb]))
    r_sat_arr = r.cartesian.without_differentials()
    v_sat_arr = r.velocity

    trial_curve = check_trial_curve(lat, lon, alt, 0, measured_curve, 0, r_sat_arr, v_sat_arr, 0, plot=True)
    diff_curve = measured_curve[:, 1] - trial_curve[:, 1]
    time = measured_curve[:, 0] - np.min(measured_curve[:, 0])
    color = curve_array[:, IDX.id] - np.min(curve_array[:, IDX.id])

    plt.figure()
    plt.scatter(time, diff_curve, c=color, label="Diff")

    fit = np.polyfit(time, diff_curve, 1)
    p = np.poly1d(fit)
    print("Fit: ", fit)

    ptime = np.linspace(np.min(time), np.max(time), 100)
    plt.plot(ptime, p(ptime), color="black", label="Fit")
    plt.figtext(.15, .12, f"Fit parameters ($Hz/s^n$): {fit}")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    find_offset_and_drift_in_data()
