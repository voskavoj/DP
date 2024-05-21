"""
    Satellite operations
"""

from astropy.time import Time
from sgp4.api import Satrec
from sgp4.conveniences import sat_epoch_datetime


class Satellite:
    """
    Satellite class
    """
    def __init__(self, constellation, number, tle_line_1, tle_line_2):
        self.name = constellation + " " + number
        try:
            self.number = int(number)
        except ValueError:
            self.number = number
        self.tle_line_1 = tle_line_1
        self.tle_line_2 = tle_line_2
        self.satrec = Satrec.twoline2rv(tle_line_1, tle_line_2)
        self.tle_time = Time(sat_epoch_datetime(self.satrec))

    def tle_age(self, time: Time):
        """
        Calculate the age of the TLE in seconds
        :param time: reference time, if None, now is used
        :return: seconds
        """
        if time is None:
            time = Time.now()
        return (time - self.tle_time).to("s").value

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def __int__(self):
        return self.number


IRIDIUM_NORAD_IDS = {
    41917: 106,
    41918: 103,
    41919: 109,
    41920: 102,
    41921: 105,
    41922: 104,
    41923: 114,
    41924: 108,
    41925: 112,
    41926: 111,
    42803: 113,
    42804: 123,
    42805: 120,
    42806: 115,
    42807: 118,
    42808: 117,
    42809: 126,
    42810: 124,
    42811: 128,
    42812: 121,
    42955: 133,
    42956: 100,
    42957: 122,
    42958: 129,
    42959: 119,
    42960: 107,
    42961: 132,
    42962: 136,
    42963: 139,
    42964: 125,
    43070: 135,
    43071: 138,
    43072: 116,
    43073: 130,
    43074: 151,
    43075: 134,
    43076: 137,
    43077: 141,
    43078: 153,
    43079: 131,
    43249: 144,
    43250: 149,
    43251: 157,
    43252: 140,
    43253: 145,
    43254: 146,
    43255: 148,
    43256: 142,
    43257: 150,
    43258: 143,
    43478: 161,
    43479: 152,
    43480: 147,
    43481: 110,
    43482: 162,
    43569: 160,
    43570: 166,
    43571: 158,
    43572: 165,
    43573: 155,
    43574: 154,
    43575: 163,
    43576: 156,
    43577: 164,
    43578: 159,
    43922: 180,
    43923: 176,
    43924: 168,
    43925: 173,
    43926: 169,
    43927: 172,
    43928: 175,
    43929: 171,
    43930: 170,
    43931: 167,
    56726: 181,
    56727: 177,
    56728: 174,
    56729: 178,
    56730: 179,
}
