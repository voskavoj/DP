"""
    Satellite operations
"""

from astropy.time import Time
from astropy.coordinates import TEME
from sgp4.api import Satrec
from sgp4.conveniences import sat_epoch_datetime

from src.satellites.predictions import predict_position_at_time


ID_NORM = {"Iridium": 10000,
           "Globalstar": 2000,
           "Orbcomm": 3000}


def generate_sattelite_id(constellation: str, sat_id: int):
    return float(ID_NORM[constellation] + sat_id)


class Satellite:
    def __init__(self, constellation, number, tle_line_1, tle_line_2):
        self.name = constellation + " " + number
        self.tle_line_1 = tle_line_1
        self.tle_line_2 = tle_line_2
        self.satrec = Satrec.twoline2rv(tle_line_1, tle_line_2)
        self.tle_time = Time(sat_epoch_datetime(self.satrec))

    def tle_age(self):
        return Time.now() - self.tle_time

    def predict_position_at_time(self, time: Time) -> TEME:
        return predict_position_at_time(self.satrec, time)
