"""
    Frame of reference transformations
"""

from astropy.time import Time
from astropy.coordinates import TEME, ITRS, CartesianDifferential, CartesianRepresentation
from astropy import units as unit


def sgp4_to_teme(pos_sgp4: tuple, vel_sgp4: tuple, time: Time):
    pos_teme = CartesianRepresentation(pos_sgp4 * unit.km)
    vel_teme = CartesianDifferential(vel_sgp4 * unit.km / unit.s)
    return TEME(pos_teme.with_differentials(vel_teme), obstime=time)


def teme_to_itrs(pos_teme: TEME) -> ITRS:
    return pos_teme.transform_to(ITRS(obstime=pos_teme.obstime))


def itrs_to_lla(pos_itrs: ITRS) -> tuple:
    location = pos_itrs.earth_location
    lat, lon, alt = location.geodetic
    return lat, lon, alt
