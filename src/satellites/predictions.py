"""
    Satellite position and observability predictions
"""

from astropy.time import Time, TimeDelta
from astropy.coordinates import TEME, EarthLocation, AltAz, ITRS, CartesianRepresentation, CartesianDifferential
from astropy import units as unit
from sgp4.api import Satrec, SGP4_ERRORS, SatrecArray
import numpy as np


def predict_satellite_positions(satrecs: list[Satrec], times: Time) -> ITRS:
    """
    Predict position in the ITRS frame of multiple satellites at multiple times using the SGP4 predictor
    For n satellites and m times, the output is an ITRS position with shape (n, m, 3)

    :param satrecs: list of satellite Satrecs
    :param times: Time object containing multiple times
    :return: ITRS positions (shape #satellites, #times, 3)
    """
    satrec_array = SatrecArray(satrecs)
    jds = np.array(times.jd1)
    jfs = np.array(times.jd2)

    errcodes, positions, velocities = satrec_array.sgp4(jds, jfs)
    pos_teme = CartesianRepresentation(x=positions[:, :, 0], y=positions[:, :, 1], z=positions[:, :, 2], unit=unit.km)
    vel_teme = CartesianDifferential(velocities[:, :, 0], velocities[:, :, 1], velocities[:, :, 2],
                                     unit=unit.km / unit.s)
    r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=times)
    pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))
    return pos_itrs


def predict_satellite_visibility(satellite: "Satellite | list", observer_position: EarthLocation, start_time: Time,
                                 duration: int = 900, step=30, azimuth_limit=None, elevation_limit=10, log=True,
                                 log_only_visible=False, native_output=False) -> list:
    """
    Predict satellite visibility

    :param satellite: Satellite object or list of Satellite objects
    :param observer_position: EarthLocation object
    :param start_time: Time object
    :param duration: duration of prediction in seconds
    :param step: step size in seconds
    :param azimuth_limit: tuple of minimum and maximum azimuth in degrees
    :param elevation_limit: minimum elevation in degrees
    :param log: print log to console
    :param log_only_visible: only log visible satellites
    :param native_output: return native times (True) or Time objects (False)

    :return: list of tuples (time, elevation, azimuth) or one tuple if only one satellite was given
    """

    # prepare array of satellites and times
    if not isinstance(satellite, list):
        satellite = [satellite]

    satrecs = SatrecArray([sat.satrec for sat in satellite])
    times = [start_time + TimeDelta(secs * unit.s) for secs in range(0, duration, step)] if duration else [start_time]
    jds = np.array([t.jd1 for t in times])
    jfs = np.array([t.jd2 for t in times])

    # do SGP4
    errcodes, positions, velocities = satrecs.sgp4(jds, jfs)

    list_of_sat_predictions = list()
    for j, sat in enumerate(satellite):
        time, errcode, pos, vel = times, errcodes[j], positions[j], velocities[j]

        if any(errcode):
            raise RuntimeError(f"SGP4 could not compute: {SGP4_ERRORS[errcode]}")

        pos_teme = CartesianRepresentation(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], unit=unit.km)
        vel_teme = CartesianDifferential(vel[:, 0], vel[:, 1], vel[:, 2], unit=unit.km / unit.s)
        r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=time)
        pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))

        topo_itrs_repr = pos_itrs.cartesian.without_differentials() - observer_position.get_itrs().cartesian
        itrs_topo = ITRS(topo_itrs_repr, obstime=time, location=observer_position)
        aa = itrs_topo.transform_to(AltAz(obstime=time, location=observer_position))
        elevs, azims = aa.alt, aa.az

        prediction_list = list()
        for i in range(len(elevs)):
            time, elevation, azimuth = times[i], elevs[i], azims[i]

            visible = True
            if elevation_limit is not None and elevation < elevation_limit * unit.deg:
                visible = False
            if azimuth_limit is not None and (azimuth_limit[0] * unit.deg <= azimuth <= azimuth_limit[1] * unit.deg):
                visible = False

            if visible:
                if native_output:
                    prediction_list.append((time.iso, elevation.value, azimuth.value))
                else:
                    prediction_list.append((time, elevation, azimuth))

        list_of_sat_predictions.append(prediction_list)

        if log:
            if len(prediction_list):
                first_time, first_ele, first_az = prediction_list[0]
                last_time, last_ele, last_az = prediction_list[-1]
                max_ele_time, max_ele_ele, max_ele_az = max(prediction_list, key=lambda x: x[1])

                print(f"{sat.name:<15}\n"
                      f"    Start: {_texttime(first_time)} @ {first_az:03.0f}\n"
                      f"    Peak:  {_texttime(max_ele_time)} @ {max_ele_az:03.0f} with elevation of {max_ele_ele:02.1f}\n"
                      f"    End:   {_texttime(last_time)} @ {last_az:03.0f} + ({(last_time - first_time).to_value('sec'):.0f} s)")
            elif not log_only_visible:
                print(f"{sat.name} will not be visible")

    if len(list_of_sat_predictions) == 1:
        return list_of_sat_predictions[0]
    else:
        return list_of_sat_predictions


def _texttime(time: Time):
    return time.value.strftime('%H:%M:%S')
