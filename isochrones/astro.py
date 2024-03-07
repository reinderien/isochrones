from datetime import datetime
from typing import NamedTuple
from zoneinfo import ZoneInfo

import numpy as np
from cartopy.crs import RotatedPole
from cartopy.feature.nightshade import Nightshade, _julian_day, _solar_position
from cartopy.geodesic import Geodesic
from shapely import Point, Polygon


# Location and time zone of the Kaaba in Mecca
KAABA_COORD = (39.826167, 21.4225)
KAABA_TIMEZONE = ZoneInfo('Asia/Riyadh')


def check_contains_night(
    local_now: datetime,
    coord: tuple[float, float],
    night: Nightshade,
) -> bool:
    geom: Polygon
    geom, = night.geometries()
    point = Point(coord)
    # This doesn't work well at all. Maybe it doesn't account for occlusion?
    return geom.contains(point) and not (6 <= local_now.hour <= 18)


def inverse_geodesic(
    coord: tuple[float, float],
    endpoint: tuple[float, float],
) -> tuple[float, float]:
    """
    Calculate the inverse geodesic parameters between the prayer location and the Kaaba. This
    uses WGS-84 and not the spherical CRS - so we use our own Geodesic instead of the spherical
    geodesic instance in Hemisphere.
    :return: Distance in metres, and departing angle in degrees counterclockwise from east.
    """
    geodesic = Geodesic()
    (distance, here_azimuth, other_azimuth), = geodesic.inverse(
        points=coord, endpoints=endpoint,
    )
    return distance, here_azimuth


class SolarPosition(NamedTuple):
    """
    hamid           nightshade
    -----           ----------
    jd              _julian_day
    d               T_UT1 (j2000)
    q               lambda_M_sun (solar longitude)
    g               M_sun (solar anomaly)
    L               lambda_ecliptic (ecliptic longitude)
    e               epsilon (ecliptic obliquity)
    D               delta_sun (declination)
    RA              alpha_sun (right ascension)
    alpha           refraction
    """
    # nightshade.py         # Hamid guide
    date: datetime          # passed into _julian_day
    T_UT1: float            # d (j2000)
    lambda_M_sun: float     # q (solar longitude)
    M_sun: float            # g (solar anomaly)
    lambda_ecliptic: float  # L (ecliptic longitude)
    epsilon: float          # e (ecliptic obliquity)
    delta_sun: float        # D (declination)
    theta_GMST: float       # Greenwich mean sidereal time (seconds)
    alpha_sun: float        # RA (right ascension) aka. lat
    pole_lat: float
    pole_lon: float         # opposite of Greenwich Hour Angle (GHA)
    central_lon: float
    rotated_pole: RotatedPole

    @classmethod
    def from_time(cls, utcnow: datetime) -> 'SolarPosition':
        """
        Adaptation of cartopy.feature.nightshade._solar_position but in rad"""

        # Centuries from J2000
        T_UT1 = (_julian_day(utcnow) - 2_451_545.0)/36_525

        # solar longitude (rad)
        lambda_M_sun = (4.894950420143297 + 628.3319872064915*T_UT1) % (2*np.pi)

        # solar anomaly (rad)
        M_sun = (6.240035938744247 + 628.3019560241842*T_UT1) % (2*np.pi)

        # ecliptic longitude (rad)
        lambda_ecliptic = (
            + lambda_M_sun
            + 0.03341723399649053 * np.sin(M_sun)
            + 3.4897235311083654e-4 * np.sin(M_sun * 2)
        )

        # obliquity of the ecliptic (epsilon in Vallado's notation)
        epsilon = 0.4090928022830742 - 2.2696610658784662e-4*T_UT1

        # declination of the sun
        delta_sun = np.arcsin(
            np.sin(epsilon) * np.sin(lambda_ecliptic)
        )

        # Greenwich mean sidereal time (seconds)
        theta_GMST = (
            + 67_310.54841
            + (876_600*3_600 + 8_640_184.812866)*T_UT1
            + 0.093104 * T_UT1**2
            - 6.2e-6 * T_UT1**3
        )

        # Convert to rad
        theta_GMST = (theta_GMST % 86_400) * 7.27220521664304e-05

        # Right ascension calculations
        numerator = (
            np.cos(epsilon) * np.sin(lambda_ecliptic)
            / np.cos(delta_sun)
        )
        denominator = np.cos(lambda_ecliptic) / np.cos(delta_sun)
        alpha_sun = np.arctan2(numerator, denominator)

        # longitude is opposite of Greenwich Hour Angle (GHA)
        # GHA == theta_GMST - alpha_sun
        lon = alpha_sun - theta_GMST
        if lon < -np.pi:
            lon += 2*np.pi

        # need longitude (opposite direction)
        lat = alpha_sun
        pole_lon = lon
        if lat > 0:
            pole_lat = lat - 0.5*np.pi
            central_lon = np.pi
        else:
            pole_lat = lat + 0.5*np.pi
            central_lon = 0

        # Skip omega0: we don't have refraction (alpha) here, and we don't care about sunrise/
        # sunset in this context

        return cls(
            date=utcnow, T_UT1=T_UT1, lambda_M_sun=lambda_M_sun, M_sun=M_sun,
            lambda_ecliptic=lambda_ecliptic, epsilon=epsilon, delta_sun=delta_sun,
            theta_GMST=theta_GMST, alpha_sun=alpha_sun,
            pole_lon=pole_lon, pole_lat=pole_lat, central_lon=central_lon,
            rotated_pole=RotatedPole(
                pole_latitude=np.rad2deg(pole_lat),
                pole_longitude=np.rad2deg(pole_lon),
                central_rotated_longitude=np.rad2deg(central_lon),
            ),
        )

    def test(self) -> None:
        comparables = _solar_position(self.date)
        assert np.allclose(
            a=np.deg2rad(comparables),
            b=(self.delta_sun, self.pole_lon),
            rtol=0, atol=1e-12,
        )

    def shadow_angle(self, shadow: float) -> float:
        # www.praytimes.org/calculation#Asr
        A = np.arccos(
            (
                np.sin(
                    np.arctan(  # arccot(1/x) = arctan(x)
                        1/(
                            shadow + np.tan(self.lambda_ecliptic - self.delta_sun)
                        )
                    )
                )
                - np.sin(self.lambda_ecliptic) * np.sin(self.delta_sun)
            ) / (
                np.cos(self.lambda_ecliptic) * np.cos(self.delta_sun)
            )
        )
        return A
