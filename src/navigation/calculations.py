import math

import numpy as np
from geopy import distance


def latlon_distance(lat1, lat2, lon1, lon2, alt1=None, alt2=None):
    """
    Returns results in m

    """
    horizontal_distance = distance.geodesic((lat1, lon1), (lat2, lon2)).m
    if alt1 is None and alt2 is None:
        return horizontal_distance
    else:
        return math.sqrt(horizontal_distance ** 2 + (alt2 - alt1) ** 2)


def m_to_deg_lat(m, lat, deg=True):
    """
    Approximate conversion from distance in km to distance in latitude
    https://en.wikipedia.org/wiki/Latitude#Meridian_distance_on_the_ellipsoid

    :param m: meters
    :param lat: current latitude (deg)
    :param deg: input is in degrees (True) or radians (False)
    :return: degrees latitude
    """
    if deg:
        lat = np.deg2rad(lat)

    m_per_deg_lat = 111132.954 - 559.822 * np.cos(2 * lat) + 1.175 * np.cos(4 * lat)
    return m / m_per_deg_lat


def m_to_deg_lon(m, lat, deg=True):
    """
    Approximate conversion from distance in km to distance in longitude
    https://en.wikipedia.org/wiki/Latitude#Meridian_distance_on_the_ellipsoid

    :param m: meters
    :param lat: current latitude (deg)
    :param deg: input is in degrees (True) or radians (False)
    :return: degrees longitude
    """
    if deg:
        lat = np.deg2rad(lat)

    m_per_deg_lon = 111132.954 * np.cos(lat)
    return m / m_per_deg_lon
