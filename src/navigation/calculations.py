import math

import numpy as np
from geopy import distance


def latlon_distance(lat1, lat2, lon1, lon2, alt1=None, alt2=None):
    """
    Calculate distance between two points in 3D space.
    :param lat1: latitude of the first point
    :param lat2: latitude of the second point
    :param lon1: longitude of the first point
    :param lon2: longitude of the second point
    :param alt1: altitude of the first point (optional)
    :param alt2: altitude of the second point (optional)
    :return: distance in meters
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


def deg_lon_to_m(lon, current_lat):
    """
    Approximate conversion from distance in longitude to distance in meters

    :param lon: longitude in degrees
    :param current_lat: current latitude in degrees
    :return: distance in meters
    """
    lat_rad = np.deg2rad(current_lat)

    m_per_deg_lon = 111132.954 * np.cos(lat_rad)
    return m_per_deg_lon * lon


def deg_lat_to_m(lat, current_lat):
    """
    Approximate conversion from distance in latitude to distance in meters

    :param lat: latitude in degrees
    :param current_lat: current latitude in degrees
    :return: distance in meters
    """
    lat_rad = np.deg2rad(current_lat)

    m_per_deg_lat = 111132.954 - 559.822 * np.cos(2 * lat_rad) + 1.175 * np.cos(4 * lat_rad)
    return m_per_deg_lat * lat


def add_gaussian_noise_and_offset(data, noise_mean=0, noise_std_dev=1, offset=0):
    """
    Add Gaussian noise and offset to the input data.

    :param data: input data
    :param noise_mean: mean of the Gaussian noise
    :param noise_std_dev: standard deviation of the Gaussian noise
    :param offset: offset to add to the data

    :return: noisy data
    """
    if noise_std_dev != 0:
        noise = np.random.normal(noise_mean, noise_std_dev, data.shape)
        data += noise

    if offset != 0:
        data += offset

    return data
