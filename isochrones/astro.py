import typing
from datetime import timedelta
from zoneinfo import ZoneInfo

import numpy as np
from cartopy.crs import RotatedPole
from cartopy.feature.nightshade import _julian_day, _solar_position

if typing.TYPE_CHECKING:
    from datetime import datetime
    from cartopy.crs import CRS, Geodetic
    from .types import (
        Coord, CoordEclipticDeg, Degree, DegArray, EclDeg, FloatArray,
        GeoDeg, Metre, Radian, Second,
    )

    TConvertInput = typing.TypeVar('TConvertInput', bound=Degree | DegArray)

    TEclipticDegree = typing.TypeVar(
        'TEclipticDegree',
        bound=DegArray | EclDeg,
    )

    class LonFromLat(typing.Protocol[TEclipticDegree]):
        def __call__(self, y_ecl: TEclipticDegree) -> TEclipticDegree:
            ...


# Location and time zone of the Kaaba in Mecca
KAABA_COORD = (39.826167, 21.4225)
KAABA_TIMEZONE = ZoneInfo('Asia/Riyadh')


def inverse_geodesic(
    geodetic: 'Geodetic', coord: 'Coord', endpoint: 'Coord',
) -> tuple[
    'Degree',  # forward azimuth, degrees counterclockwise from east
    'Degree',  # back azimuth, degrees counterclockwise from east
    'Metre',   # distance, m
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
    'Degree',  # forward azimuth, degrees counterclockwise from east
    'Degree',  # back azimuth, degrees counterclockwise from east
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


def unwrap(p: 'Radian') -> 'Radian':
    """
    Vaguely similar to numpy.unwrap, but producing always-positive angles, and
    operating on a scalar.
    """
    return p % (2*np.pi)


def shadow_angle(y_ecl: 'DegArray', shadow: float) -> 'DegArray':
    """
    Based on https://radhifadlillah.com/blog/2020-09-06-calculating-prayer-times/
    Return the angular difference in rad between solar noon and the 'asr'
    """
    # The reference implementations rely on delta_sun (declination). Since
    # we defer declination calculation to a transformation after this
    # function, it's already accounted for, and any declination terms are
    # replaced with sin(0)=0 and cos(0)=1.
    y_rad = np.deg2rad(y_ecl)
    arg = np.sin(
        np.arctan(  # arccot(1/x) = arctan(x)
            1/(
                shadow + np.tan(
                    np.abs(y_rad)
                )
            )
        )
    ) / np.cos(y_rad)

    # Let the NaNs through, but don't complain about them.
    # arg = np.clip(arg, a_min=-1, a_max=1)
    is_valid = (arg >= -1) & (arg <= +1)
    A: FloatArray = np.full_like(a=arg, fill_value=np.nan)
    A[is_valid] = np.arccos(arg[is_valid])

    return np.rad2deg(A)


def convert_crs(
    source: 'CRS', dest: 'CRS',
    x: 'TConvertInput',
    y: 'TConvertInput',
) -> 'DegArray':
    if isinstance(x, float) and isinstance(y, float):
        x_array: 'DegArray' = np.array((x,))
        y_array: 'DegArray' = np.array((y,))
    else:
        x_array = x
        y_array = y
    xyz: DegArray = dest.transform_points(
        x=x_array, y=y_array, src_crs=source,
    ).T
    return xyz[:-1, ...]


def refraction_to_ecliptic_longitude(
    y_ecl: 'typing.Union[EclDeg, FloatArray]',
    refraction_deg: 'EclDeg',
    pm: bool,
) -> 'typing.Union[EclDeg, FloatArray]':
    # Adaptation of Nightshade(), but with a coordinate system fixup, and only populating one
    # side of daytime
    refraction = np.deg2rad(refraction_deg)

    # Solve the generalized equation for omega0, which is the
    # angle of sunrise/sunset from solar noon
    # We need to clip the input to arccos to [-1, 1] due to floating
    # point precision and arccos creating nans for values outside
    # of the domain
    arccos_tmp = np.clip(
        a=np.sin(refraction) / np.cos(np.deg2rad(y_ecl)),
        a_min=-1, a_max=1,
    )
    omega0 = np.arccos(arccos_tmp)

    # Fill the longitude values from the offset for midnight.
    x = omega0 - np.pi
    if not pm:
        x = -x

    return np.rad2deg(x)


DEFAULT_Y_ECLIPTIC: 'DegArray' = np.linspace(
    start=-90, stop=+90, num=91, dtype=np.float64,
)


class SolarPosition(typing.NamedTuple):
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
    # nightshade.py            # Hamid guide
    date: 'datetime'           # passed into _julian_day
    T_UT1: float               # d (j2000)
    lambda_M_sun: 'Radian'     # q (solar longitude)
    M_sun: 'Radian'            # g (solar anomaly)
    lambda_ecliptic: 'Radian'  # L (ecliptic longitude)
    epsilon: 'Radian'          # e (ecliptic obliquity)
    delta_sun: 'Radian'        # D (declination)
    theta_GMST: 'Second'       # Greenwich mean sidereal time (seconds)
    alpha_sun: 'Radian'        # RA (right ascension) aka. lat
    pole_lat: 'Radian'
    pole_lon: 'Radian'         # opposite of Greenwich Hour Angle (GHA)
    central_lon: 'Radian'
    rotated_pole: RotatedPole

    @classmethod
    def from_time(cls, utcnow: 'datetime') -> 'SolarPosition':
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

    def ecliptic_to_geographic(
        self,
        globe_crs: 'CRS',
        x_ecl,
        y_ecl,
    ) -> 'FloatArray':
        return convert_crs(
            source=globe_crs, dest=self.rotated_pole, x=x_ecl, y=y_ecl,
        )

    def geographic_to_ecliptic(
        self,
        globe_crs: 'CRS',
        x_geo,
        y_geo,
    ) -> 'FloatArray':
        return convert_crs(
            source=self.rotated_pole, dest=globe_crs, x=x_geo, y=y_geo,
        )

    def ecliptic_to_geolon_time(
        self,
        globe_crs: 'CRS',
        xy: 'CoordEclipticDeg',
    ) -> tuple['GeoDeg', timedelta]:
        x_ecl, y_ecl = xy

        # This is solar time. We need coordinated time eventually, which will add the
        # current local time and the scaled difference between ecliptic longitudes.
        time = timedelta(days=x_ecl/360 + 0.5)

        x_geo, y_geo = self.ecliptic_to_geographic(
            globe_crs=globe_crs, x_ecl=x_ecl, y_ecl=y_ecl,
        )
        return x_geo, time

    def noon_angle_to_geolon_time(
        self,
        globe_crs: 'CRS',
        make_lon: 'LonFromLat[EclDeg]',
        y_ecl: 'EclDeg',
    ) -> tuple['GeoDeg', timedelta]:
        x_ecl = make_lon(y_ecl=y_ecl)
        return self.ecliptic_to_geolon_time(
            globe_crs=globe_crs, xy=(x_ecl, y_ecl),
        )

    def noon_angle_to_isochrone(
        self,
        globe_crs: 'CRS',
        make_lon: 'LonFromLat[FloatArray]',
    ) -> 'FloatArray':
        """
        :param globe_crs: The coordinate reference system of the globe, used when translating to the
                          night-rotated coordinate system. Typically Geodetic.
        :return: A 2*n array of x, y coordinates in degrees; the isochrone curve.
        """
        x_ecl = make_lon(y_ecl=DEFAULT_Y_ECLIPTIC)
        return self.ecliptic_to_geographic(
            globe_crs=globe_crs, x_ecl=x_ecl, y_ecl=DEFAULT_Y_ECLIPTIC,
        )

    def refraction_to_geolon_time(
        self,
        globe_crs: 'CRS',
        y_ecl: 'EclDeg',
        refraction: 'Degree',
        pm: bool,
    ) -> tuple['GeoDeg', timedelta]:
        x_ecl = refraction_to_ecliptic_longitude(
            refraction_deg=refraction, pm=pm, y_ecl=y_ecl,
        )
        return self.ecliptic_to_geolon_time(
            globe_crs=globe_crs, xy=(x_ecl, y_ecl),
        )

    def refraction_to_isochrone(
        self,
        globe_crs: 'CRS',
        refraction: 'Degree',
        pm: bool,
    ) -> 'FloatArray':
        x_ecl = refraction_to_ecliptic_longitude(
            refraction_deg=refraction, pm=pm, y_ecl=DEFAULT_Y_ECLIPTIC,
        )
        return self.ecliptic_to_geographic(
            globe_crs=globe_crs, x_ecl=x_ecl, y_ecl=DEFAULT_Y_ECLIPTIC,
        )
