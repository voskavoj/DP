"""
    Satellite position and observability predictions
"""

from astropy.time import Time
from astropy.coordinates import TEME, EarthLocation, AltAz, ITRS
from astropy import units as unit
from sgp4.api import Satrec, SGP4_ERRORS

from src.satellites.coord_transforms import teme_to_itrs, itrs_to_lla, sgp4_to_teme


def predict_position_at_time(satrec: Satrec, time: Time) -> TEME:
    errcode, pos, vel = satrec.sgp4(time.jd1, time.jd2)

    if errcode:
        raise RuntimeError(f"SGP4 could not compute: {SGP4_ERRORS[errcode]}")

    return sgp4_to_teme(pos, vel, time)


def predict_satellite_visibility(satellite: 'Satellite', observer_position: EarthLocation, start_time: Time,
                                 duration: int = 900, step=30, azimuth_limit=None, elevation_limit=10, log=True,
                                 log_only_visible=False):
    prediction_list = list()
    for i in range(duration // step):
        time = start_time + (i * step * unit.s)
        pos_itrs = teme_to_itrs(predict_position_at_time(satellite.satrec, time))

        topo_itrs_repr = pos_itrs.cartesian.without_differentials() - observer_position.get_itrs(time).cartesian
        itrs_topo = ITRS(topo_itrs_repr, obstime=time, location=observer_position)
        aa = itrs_topo.transform_to(AltAz(obstime=time, location=observer_position))
        elevation, azimuth = aa.alt, aa.az

        visible = True
        if elevation_limit is not None and elevation < elevation_limit * unit.deg:
            visible = False
        if azimuth_limit is not None and (azimuth_limit[0] * unit.deg <= azimuth <= azimuth_limit[1] * unit.deg):
            visible = False

        if visible:
            prediction_list.append((time, elevation, azimuth))

    if log:
        if len(prediction_list):
            first_time, first_ele, first_az = prediction_list[0]
            last_time, last_ele, last_az = prediction_list[-1]
            max_ele_time, max_ele_ele, max_ele_az = max(prediction_list, key=lambda x: x[1])

            print(f"{satellite.name:<15}\n"
                  f"    Start: {_texttime(first_time)} @ {first_az:03.0f}\n"
                  f"    Peak:  {_texttime(max_ele_time)} @ {max_ele_az:03.0f} with elevation of {max_ele_ele:02.1f}\n"
                  f"    End:   {_texttime(last_time)} @ {last_az:03.0f} + ({(last_time - first_time).to_value('sec'):.0f} s)")
        elif not log_only_visible:
            print(f"{satellite.name} will not be visible")

    if len(prediction_list):
        return prediction_list
    else:
        return False


def _texttime(time: Time):
    return time.value.strftime('%H:%M:%S')
