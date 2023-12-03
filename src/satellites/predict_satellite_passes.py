"""
    App
    Predicts close satellite passes
"""
from astropy.coordinates import EarthLocation
from astropy.time import Time

from src.satellites.download_tle import download_tles, unpack

# Constants
from src.satellites.predictions import predict_satellite_visibility

PREDICTION_MINUTES = 15

#                   lon (east)   lat (north) alt (m)
LOCATIONS = {"FEL": (14.3923333, 50.1030019, 225),
             }

CONSTELLATIONS = ("Iridium", )


# App
observer_position = EarthLocation.from_geodetic(*LOCATIONS["FEL"])
satellites = unpack(download_tles(CONSTELLATIONS))
time = Time.now()

for sat in satellites:
    t = Time.now()
    predict_satellite_visibility(sat, observer_position, time,
                                 duration=PREDICTION_MINUTES * 60, step=120,
                                 log=True, log_only_visible=True)
    print((Time.now() - t).to_value("sec"))



