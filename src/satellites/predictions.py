"""
    Satellite position and observability predictions
"""

from astropy.time import Time, TimeDelta
from astropy.coordinates import TEME, EarthLocation, AltAz, ITRS, CartesianRepresentation, CartesianDifferential
from astropy import units as unit
from sgp4.api import Satrec, SGP4_ERRORS, SatrecArray
import numpy as np

from src.satellites.coord_transforms import teme_to_itrs, itrs_to_lla, sgp4_to_teme


def predict_position_at_time(satrec: Satrec, time: Time) -> TEME:
    errcode, pos, vel = satrec.sgp4(time.jd1, time.jd2)

    if errcode:
        raise RuntimeError(f"SGP4 could not compute: {SGP4_ERRORS[errcode]}")

    return sgp4_to_teme(pos, vel, time)


def predict_satellite_visibility(satellite: "Satellite | list", observer_position: EarthLocation, start_time: Time,
                                 duration: int = 900, step=30, azimuth_limit=None, elevation_limit=10, log=True,
                                 log_only_visible=False) -> list:

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
