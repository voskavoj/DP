"""
    App
    Predicts close satellite passes
"""
from astropy.coordinates import EarthLocation
from astropy.time import Time
from time import sleep

from src.satellites.download_tle import download_tles, unpack

from src.config.setup import *
from src.satellites.predictions import predict_satellite_visibility, _texttime
from src.config.locations import LOCATIONS  # (lon (east), lat (north), alt (m))

PREDICTION_MINUTES = 0
TIME = "2024-01-07 18:17:40"  # UTC


# App
observer_position = EarthLocation.from_geodetic(*LOCATIONS[LOCATION])
satellites = unpack(download_tles(CONSTELLATIONS))
time = Time(TIME, format="iso") if TIME else Time.now()

print(f"Prediction {PREDICTION_MINUTES} minutes ahead")

predict_satellite_visibility(satellites, observer_position, time,
                             duration=PREDICTION_MINUTES * 60, step=5,
                             elevation_limit=10,
                             log=True, log_only_visible=True)

while True:
    time = Time.now()
    pred_list = list()
    prediction = predict_satellite_visibility(satellites, observer_position, time,
                                              duration=0, step=5,
                                              elevation_limit=10,
                                              log=False)

    for i, sat in enumerate(satellites):
        pred = prediction[i]
        if pred:
            _, ele, az = pred[0]
            pred_list.append((sat.name, ele, az))

    print(f"{_texttime(time)}: {' ||| '.join(f'{name:<10} @ {az:03.0f} @ {ele:02.1f}' for name, ele, az in pred_list)}")
    sleep(2)
