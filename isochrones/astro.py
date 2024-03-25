import typing
from datetime import datetime
from typing import NamedTuple
from zoneinfo import ZoneInfo

import numpy as np
from cartopy.crs import RotatedPole
from cartopy.feature.nightshade import _julian_day, _solar_position

if typing.TYPE_CHECKING:
    from cartopy.crs import CRS, Geodetic
    from .types import Coord, FloatArray

    class LonFromLat(typing.Protocol):
        def __call__(self, y: 'FloatArray') -> 'FloatArray':
            ...


# Location and time zone of the Kaaba in Mecca
KAABA_COORD = (39.826167, 21.4225)
KAABA_TIMEZONE = ZoneInfo('Asia/Riyadh')


def inverse_geodesic(
    geodetic: 'Geodetic', coord: 'Coord', endpoint: 'Coord',
) -> tuple[
    float,  # forward azimuth, degrees counterclockwise from east
    float,  # back azimuth, degrees counterclockwise from east
    float,  # distance, m
]:
    """
    Calculate the inverse geodesic parameters between the prayer location and the Kaaba.
    Validate with e.g. https://geographiclib.sourceforge.io/cgi-bin/GeodSolve
    """
    geodesic = geodetic.get_geod()

    # unpack-repack to help out the type system
    forward, back, distance = geodesic.inv(
        lons1=coord[0], lons2=endpoint[0],
        lats1=coord[1], lats2=endpoint[1],
    )
    return forward % 360, back % 360, distance


def inverse_geodesic_hack(
    sphere: 'Geodetic', coord: 'Coord', endpoint: 'Coord',
) -> tuple[
    float,  # forward azimuth, degrees counterclockwise from east
    float,  # back azimuth, degrees counterclockwise from east
]:
    # This is horrible, but it's the only way I can get the geodesic to match the
    # departing angle on the orthographic plot
    result = sphere.get_geod().inv_intermediate(
        lon1=coord[0], lat1=coord[1],
        lon2=endpoint[0], lat2=endpoint[1],
        npts=int(endpoint[0] - coord[0]) // 10,
        return_back_azimuth=True,
    )
    lons = np.array((result.lons[0], result.lons[-1])) - (coord[0], endpoint[0])
    lats = np.array((result.lats[0], result.lats[-1])) - (coord[1], endpoint[1])
    forward, back = np.rad2deg(np.arctan2(lats, lons)) % 360
    return forward, back


def unwrap(p: float) -> float:
    """
    Vaguely similar to numpy.unwrap, but producing always-positive angles, and
    operating on a scalar.
    """
    return p % (2*np.pi)


def shadow_angle(y: 'FloatArray', shadow: float) -> 'FloatArray':
    """
    Based on https://radhifadlillah.com/blog/2020-09-06-calculating-prayer-times/
    Return the angular difference in rad between solar noon and the 'asr'
    """
    # The reference implementations rely on delta_sun (declination). Since
    # we defer declination calculation to a transformation after this
    # function, it's already accounted for, and any declination terms are
    # replaced with sin(0)=0 and cos(0)=1.
    arg = np.sin(
        np.arctan(  # arccot(1/x) = arctan(x)
            1/(
                shadow + np.tan(
                    np.abs(y)
                )
            )
        )
    ) / np.cos(y)

    # Let the NaNs through, but don't complain about them.
    # arg = np.clip(arg, a_min=-1, a_max=1)
    is_valid = (arg >= -1) & (arg <= +1)
    A: FloatArray = np.full_like(a=arg, fill_value=np.nan)
    A[is_valid] = np.arccos(arg[is_valid])

    return A


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
        Adaptation of cartopy.feature.nightshade._solar_position but in rad
        """

        # Centuries from J2000
        T_UT1 = (_julian_day(utcnow) - 2_451_545.)/36_525

        # solar longitude (rad)
        lambda_M_sun = unwrap(4.894950420143297 + 628.3319872064915*T_UT1)

        # solar anomaly (rad)
        M_sun = unwrap(6.240035938744247 + 628.3019560241842*T_UT1)

        # ecliptic longitude (rad)
        lambda_ecliptic = unwrap(
            + lambda_M_sun
            + 0.03341723399649053 * np.sin(M_sun)
            + 3.4897235311083654e-4 * np.sin(M_sun * 2)
        )

        # obliquity of the ecliptic (epsilon in Vallado's notation)
        epsilon = unwrap(0.4090928022830742 - 2.2696610658784662e-4*T_UT1)

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
        theta_GMST = unwrap((theta_GMST % 86_400) * 7.27220521664304e-05)

        # Right ascension calculations
        numerator = np.cos(epsilon) * np.sin(lambda_ecliptic)
        denominator = np.cos(lambda_ecliptic)
        alpha_sun = np.arctan2(numerator, denominator)

        # longitude is opposite of Greenwich Hour Angle (GHA)
        # GHA == theta_GMST - alpha_sun
        lon = alpha_sun - theta_GMST
        if lon < -np.pi:
            lon += 2*np.pi

        # need longitude (opposite direction)
        lat = delta_sun
        pole_lon = lon

        # Original Nightshade() has a conditional here that changes the angle addend to lat.
        # That produces incorrect behaviour - an inverted coordinate system.
        pole_lat = lat + 0.5*np.pi
        central_lon = 0

        # Postpone calculation of omega0 until isochrone routine: we don't have refraction
        # (alpha) here, and we don't care about sunrise/sunset in this context

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

    def isochrone_from_noon_angle(
        self, globe_crs: 'CRS', make_lon: 'LonFromLat',
    ) -> 'FloatArray':
        """
        :param globe_crs: The coordinate reference system of the globe, used when translating to the
                          night-rotated coordinate system. Typically Geodetic.
        :return: A 2*n array of x, y coordinates in degrees; the isochrone curve.
        """
        y = np.linspace(start=-np.pi/2, stop=+np.pi/2, num=91)
        x = make_lon(y=y)
        xyz: FloatArray = globe_crs.transform_points(
            x=np.rad2deg(x) + 180,
            y=np.rad2deg(y),
            src_crs=self.rotated_pole,
        ).T
        return xyz[:-1]

    def isochrone_from_refraction(
        self, globe_crs: 'CRS', angle: float, pm: bool,
    ) -> 'FloatArray':
        # Adaptation of Nightshade(), but with a coordinate system fixup, and only populating one
        # side of daytime
        refraction = np.deg2rad(angle)
        npts = 91

        # Fill latitudes up
        y = np.linspace(-0.5*np.pi - refraction, 0.5*np.pi + refraction, npts)

        # Solve the generalized equation for omega0, which is the
        # angle of sunrise/sunset from solar noon
        # We need to clip the input to arccos to [-1, 1] due to floating
        # point precision and arccos creating nans for values outside
        # of the domain
        arccos_tmp = np.clip(
            a=np.sin(refraction) / np.cos(y),
            a_min=-1, a_max=1,
        )
        omega0 = np.arccos(arccos_tmp)

        # Fill the longitude values from the offset for midnight.
        x = omega0 - np.pi
        if not pm:
            x = -x

        xyz: FloatArray = globe_crs.transform_points(
            x=np.rad2deg(x),
            y=np.rad2deg(y),
            src_crs=self.rotated_pole,
        ).T
        return xyz[:-1]
