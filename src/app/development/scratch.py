from astropy.coordinates import CartesianRepresentation, CartesianDifferential, TEME, ITRS
from astropy.time import Time
from astropy import units as unit

from src.config.setup import CONSTELLATIONS
from src.satellites.download_tle import download_tles
from src.navigation.calculations import latlon_distance


def verify_sat_pos_vs_orbitron():
    time = Time("2024-04-05 15:00:00", format="iso", scale="utc")
    satellites = download_tles(CONSTELLATIONS)
    sat = satellites["Iridium"]["103"].satrec

    jd, jf = time.jd1, time.jd2
    err, pos, vel = sat.sgp4(jd, jf)
    pos_teme = CartesianRepresentation(pos, unit=unit.km)
    vel_teme = CartesianDifferential(vel, unit=unit.km / unit.s)
    r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=time)
    pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))

    geo_pos = pos_itrs.earth_location.geodetic

    print(geo_pos)


verify_sat_pos_vs_orbitron()
print(latlon_distance(90 - 167.61238366, 90 - 167.6119, 77.76230298, 77.7623, 789.14026782, 789.142))
