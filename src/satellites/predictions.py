"""
    Satellite position and observability predictions
"""

from astropy.time import Time, TimeDelta
from astropy.coordinates import TEME, EarthLocation, AltAz, ITRS, CartesianRepresentation, CartesianDifferential
from astropy import units as unit
from sgp4.api import Satrec, SGP4_ERRORS, SatrecArray
import numpy as np

from src.navigation.calculations import latlon_distance

C = 299792458


def predict_satellite_positions(satrecs: list[Satrec], times: Time) -> ITRS:
    """
    Predict position in the ITRS frame of multiple satellites at multiple times using the SGP4 predictor
    For n satellites and m times, the output is an ITRS position with shape (n, m, 3)

    :param satrecs: list of satellite Satrecs
    :param times: Time object containing multiple times
    :return: ITRS positions (shape #satellites, #times, 3)
    """
    satrec_array = SatrecArray(satrecs)
    jds = np.array(times.jd1)
    jfs = np.array(times.jd2)

    errcodes, positions, velocities = satrec_array.sgp4(jds, jfs)
    pos_teme = CartesianRepresentation(x=positions[:, :, 0], y=positions[:, :, 1], z=positions[:, :, 2], unit=unit.km)
    vel_teme = CartesianDifferential(velocities[:, :, 0], velocities[:, :, 1], velocities[:, :, 2],
                                     unit=unit.km / unit.s)
    r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=times)
    pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))
    return pos_itrs


def predict_satellite_doppler_shift(satellite: "Satellite | list", lat: float, lon: float, alt: float, base_freq: int,
                                    start_time: Time, duration: int, step: int | float, elevation_limit=0) -> np.array:
    """
    Predict doppler shift of multiple satellites
    Returns array [rel_time, shifted frequency, abs_doppler_shift, rel_doppler_shift]
    If elevation_limit is set, invalid values will be [0, 0, 0, 0]

    :param satellite: Satellite or list of Satellites
    :param lat: Observer latitude
    :param lon: Observer longitude
    :param alt: Observer altitude
    :param base_freq: base frequency in Hz
    :param start_time: start time as Time
    :param duration: duration of prediction in seconds
    :param step: step size in seconds (supports floats)

    :return: array of predictions of shape #sat x #times x 4
    """

    c = C * unit.m / unit.s
    base_freq = base_freq * unit.Hz

    # prepare array of satellites and times
    if not isinstance(satellite, list):
        satellite = [satellite]

    satrecs = SatrecArray([sat.satrec for sat in satellite])
    times = [start_time + TimeDelta(secs * unit.s)
             for secs in [round(i * step, 5) for i in range(int(duration/step))]] if duration else [start_time]
    jds = np.array([t.jd1 for t in times])
    jfs = np.array([t.jd2 for t in times])

    observer_position = EarthLocation.from_geodetic([lon] * len(times), [lat] * len(times), [alt] * len(times))

    # do SGP4
    errcodes, positions, velocities = satrecs.sgp4(jds, jfs)

    pred_array = np.empty((len(satellite), len(times) - 1, 4))
    for j, sat in enumerate(satellite):
        time, errcode, pos, vel = Time(times), errcodes[j], positions[j], velocities[j]

        if any(errcode):
            raise RuntimeError(f"SGP4 could not compute: {SGP4_ERRORS[errcode]}")

        pos_teme = CartesianRepresentation(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], unit=unit.km)
        vel_teme = CartesianDifferential(vel[:, 0], vel[:, 1], vel[:, 2], unit=unit.km / unit.s)
        r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=time)
        pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))

        if elevation_limit is not None:
            topo_itrs_repr = pos_itrs.cartesian.without_differentials() - observer_position.get_itrs(
                obstime=time).cartesian
            itrs_topo = ITRS(topo_itrs_repr, obstime=time, location=observer_position)
            aa = itrs_topo.transform_to(AltAz(obstime=time, location=observer_position))
            elevs = aa.alt

        user_pos_itrs = observer_position.get_itrs(time)
        distance = pos_itrs.cartesian.without_differentials() - user_pos_itrs.cartesian.without_differentials()
        distance = distance.norm()

        for i in range(len(distance) - 1):
            if elevation_limit is not None and elevs[i] < elevation_limit:
                continue

            delta_distance = (distance[i] - distance[i + 1]).to(unit.m)
            delta_time = (time[i + 1] - time[i]).sec * unit.s
            rel_velocity = delta_distance / delta_time
            rel_doppler_shift = c / (c - rel_velocity)
            abs_doppler_shift = (rel_doppler_shift * base_freq) - base_freq
            shifted_freq = base_freq + abs_doppler_shift

            # output array: rel_time, shifted frequency, abs_doppler_shift, rel_doppler_shift
            pred_array[j, i] = [(time[i] - time[0]).sec,
                                shifted_freq.value, abs_doppler_shift.value, rel_doppler_shift.value]

    return pred_array


def predict_satellite_visibility(satellite: "Satellite | list", observer_position: EarthLocation, start_time: Time,
                                 duration: int = 900, step=30, azimuth_limit=None, elevation_limit=10, log=True,
                                 log_only_visible=False, native_output=False) -> list:
    """
    Predict satellite visibility

    :param satellite: Satellite object or list of Satellite objects
    :param observer_position: EarthLocation object
    :param start_time: Time object
    :param duration: duration of prediction in seconds
    :param step: step size in seconds
    :param azimuth_limit: tuple of minimum and maximum azimuth in degrees
    :param elevation_limit: minimum elevation in degrees
    :param log: print log to console
    :param log_only_visible: only log visible satellites
    :param native_output: return native times (True) or Time objects (False)

    :return: list of tuples (time, elevation, azimuth) or one tuple if only one satellite was given
    """

    # prepare array of satellites and times
    if not isinstance(satellite, list):
        satellite = [satellite]

    satrecs = SatrecArray([sat.satrec for sat in satellite])
    times = [start_time + TimeDelta(secs * unit.s) for secs in range(0, duration, step)] if duration else [start_time]
    jds = np.array([t.jd1 for t in times])
    jfs = np.array([t.jd2 for t in times])

    # do SGP4
    errcodes, positions, velocities = satrecs.sgp4(jds, jfs)

    list_of_sat_predictions = list()
    for j, sat in enumerate(satellite):
        time, errcode, pos, vel = times, errcodes[j], positions[j], velocities[j]

        if any(errcode):
            raise RuntimeError(f"SGP4 could not compute: {SGP4_ERRORS[errcode]}")

        pos_teme = CartesianRepresentation(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], unit=unit.km)
        vel_teme = CartesianDifferential(vel[:, 0], vel[:, 1], vel[:, 2], unit=unit.km / unit.s)
        r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=time)
        pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))

        topo_itrs_repr = pos_itrs.cartesian.without_differentials() - observer_position.get_itrs().cartesian
        itrs_topo = ITRS(topo_itrs_repr, obstime=time, location=observer_position)
        aa = itrs_topo.transform_to(AltAz(obstime=time, location=observer_position))
        elevs, azims = aa.alt, aa.az

        prediction_list = list()
        for i in range(len(elevs)):
            time, elevation, azimuth = times[i], elevs[i], azims[i]

            visible = True
            if elevation_limit is not None and elevation < elevation_limit * unit.deg:
                visible = False
            if azimuth_limit is not None and (azimuth_limit[0] * unit.deg <= azimuth <= azimuth_limit[1] * unit.deg):
                visible = False

            if visible:
                if native_output:
                    prediction_list.append((time.iso, elevation.value, azimuth.value))
                else:
                    prediction_list.append((time, elevation, azimuth))

        list_of_sat_predictions.append(prediction_list)

        if log:
            if len(prediction_list):
                first_time, first_ele, first_az = prediction_list[0]
                last_time, last_ele, last_az = prediction_list[-1]
                max_ele_time, max_ele_ele, max_ele_az = max(prediction_list, key=lambda x: x[1])

                print(f"{sat.name:<15}\n"
                      f"    Start: {_texttime(first_time)} @ {first_az:03.0f}\n"
                      f"    Peak:  {_texttime(max_ele_time)} @ {max_ele_az:03.0f} with elevation of {max_ele_ele:02.1f}\n"
                      f"    End:   {_texttime(last_time)} @ {last_az:03.0f} + ({(last_time - first_time).to_value('sec'):.0f} s)")
            elif not log_only_visible:
                print(f"{sat.name} will not be visible")

    if len(list_of_sat_predictions) == 1:
        return list_of_sat_predictions[0]
    else:
        return list_of_sat_predictions


def find_closest_tle_id_to_ira_id_by_position(time, satrecs, sat_list, iri_lat, iri_lon, iri_alt):
    """
    Find the closest satellite to the Iridium satellite by position
    :param time: reference time
    :param satrecs: satrecarray
    :param sat_list: list of TLEs
    :param iri_lat: latitude from IRA frame
    :param iri_lon: longitude from IRA frame
    :param iri_alt: altitude from IRA frame
    :return: closest satellite ID, distance, SGP4 latitude, longitude, altitude
    """

    jds = np.array([time.jd1])
    jfs = np.array([time.jd2])
    errcodes, positions, velocities = satrecs.sgp4(jds, jfs)
    pos_teme = CartesianRepresentation(x=positions[:, 0, 0], y=positions[:, 0, 1], z=positions[:, 0, 2], unit=unit.km)
    vel_teme = CartesianDifferential(velocities[:, 0, 0], velocities[:, 0, 1], velocities[:, 0, 2],
                                     unit=unit.km / unit.s)
    r_teme = TEME(pos_teme.with_differentials(vel_teme), obstime=time)
    pos_itrs = r_teme.transform_to(ITRS(obstime=r_teme.obstime))
    geo_pos = pos_itrs.earth_location.geodetic

    distances = [latlon_distance(iri_lat, lat.value, iri_lon, lon.value, iri_alt*1000, alt.value*1000) for lat, lon, alt in
                 zip(geo_pos.lat, geo_pos.lon, geo_pos.height.to("km"))]
    min_idx = np.argmin(np.array(distances))
    closest_sat_id = sat_list[min_idx].number

    dist = distances[min_idx]
    tle_lat, tle_lon, tle_alt = geo_pos.lat[min_idx].value, geo_pos.lon[min_idx].value, geo_pos.height[min_idx].value

    return closest_sat_id, dist, tle_lat, tle_lon, tle_alt


def _texttime(time: Time):
    """
        Deprecated
    """
    return time
