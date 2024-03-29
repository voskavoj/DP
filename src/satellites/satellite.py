"""
    Satellite operations
"""

from astropy.time import Time
from sgp4.api import Satrec
from sgp4.conveniences import sat_epoch_datetime


class Satellite:
    def __init__(self, constellation, number, tle_line_1, tle_line_2):
        self.name = constellation + " " + number
        self.number = int(number)
        self.tle_line_1 = tle_line_1
        self.tle_line_2 = tle_line_2
        self.satrec = Satrec.twoline2rv(tle_line_1, tle_line_2)
        self.tle_time = Time(sat_epoch_datetime(self.satrec))

    def tle_age(self):
        return Time.now() - self.tle_time

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def __int__(self):
        return self.number
