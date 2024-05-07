import numpy as np
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, TEME, ITRS
from astropy.time import Time
from astropy import units as unit

from src.config.setup import CONSTELLATIONS
from src.radio.iridium_channels import IRI_ID_TO_TLE_ID_MAP
from src.satellites.download_tle import download_tles
from src.navigation.calculations import latlon_distance, deg_lat_to_m, deg_lon_to_m


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


# verify_sat_pos_vs_orbitron()
# print(latlon_distance(90 - 167.61238366, 90 - 167.6119, 77.76230298, 77.7623, 789.14026782, 789.142))

print(f"{deg_lat_to_m(1,  0)/1000:.3f}, {deg_lon_to_m(1,  0)/1000:.3f}")
print(f"{deg_lat_to_m(1, 15)/1000:.3f}, {deg_lon_to_m(1, 15)/1000:.3f}")
print(f"{deg_lat_to_m(1, 30)/1000:.3f}, {deg_lon_to_m(1, 30)/1000:.3f}")
print(f"{deg_lat_to_m(1, 45)/1000:.3f}, {deg_lon_to_m(1, 45)/1000:.3f}")
print(f"{deg_lat_to_m(1, 60)/1000:.3f}, {deg_lon_to_m(1, 60)/1000:.3f}")
print(f"{deg_lat_to_m(1, 75)/1000:.3f}, {deg_lon_to_m(1, 75)/1000:.3f}")
print(f"{deg_lat_to_m(1, 90)/1000:.3f}, {deg_lon_to_m(1, 90)/1000:.3f}")


print()
print(f"{deg_lat_to_m(1,  50.1)/1000:.3f}, {deg_lon_to_m(1,  50.1)/1000:.3f}")
print(f"{deg_lat_to_m(0.01,  50.1)/1000:.3f}, {deg_lon_to_m(0.01,  50.1)/1000:.3f}")


id_table = np.loadtxt("Data\\sat08\\id.csv", dtype="int", delimiter=";")
# (norad, frame)
print(id_table.shape)
for idx_iri in range(id_table.shape[1]):
    i1, c1, i2, c2 = np.argmax(id_table[:, idx_iri]), np.max(id_table[:, idx_iri]), np.argsort(id_table[:, idx_iri])[-2], np.sort(id_table[:, idx_iri])[-2]
    missed_ratio = ((np.sum(id_table[:, idx_iri]) - c1) / c1) if c1 else 0
    outstr = ""
    if c1:
        outstr = f"{idx_iri:03d} & {i1:03d} & {c1: 5d} &"
        if c2:
            outstr += f" {i2:03d} & {c2: 5d} &"
        else:
            outstr += "     &       &"

        outstr += f" {100 - missed_ratio*100:.0f}\%     \\\\"

        print(outstr)

print()
for idx_iri in range(id_table.shape[1]):
    i1, c1, i2, c2 = np.argmax(id_table[:, idx_iri]), np.max(id_table[:, idx_iri]), np.argsort(id_table[:, idx_iri])[-2], np.sort(id_table[:, idx_iri])[-2]
    new_id = IRI_ID_TO_TLE_ID_MAP.get(idx_iri, None)

    if new_id is not None:
        # outstr = ""
        # if c1:
        outstr = f"{idx_iri:03d} & {i1:03d} & {new_id:03d} & {new_id:03d} \\\\"
        print(outstr)


import astropy
print(astropy.__citation__)