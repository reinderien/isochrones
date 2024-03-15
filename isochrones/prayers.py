from array import array
from dataclasses import dataclass
from datetime import datetime

import numpy as np
from cartopy.crs import CRS
from cartopy.feature.nightshade import Nightshade

from .astro import SolarPosition


def isochrone_from_noon_angle(globe_crs: CRS, utcnow: datetime, angle: float) -> np.ndarray:
    """
    :param globe_crs: The coordinate reference system of the globe, used when translating to the
                      night-rotated coordinate system. Typically Geodetic.
    :param utcnow: Timezone-aware datetime used to locate the sun
    :param angle: Degrees after noon
    :return: A 2*n array of x, y coordinates in degrees; the isochrone curve.
    """
    sun = SolarPosition.from_time(utcnow=utcnow)
    sun.test()
    y = np.linspace(start=-90, stop=+90, num=91)
    x = np.full_like(a=y, fill_value=180 + angle)
    xyz = globe_crs.transform_points(
        x=x, y=y, src_crs=sun.rotated_pole,
    ).T
    return xyz[:-1]


@dataclass(frozen=True)
class Prayer:
    """One salah prayer definition"""
    name: str     # prayer name in romanized Arabic
    colour: str   # time-invariant graphing colour

    def isochrone(self, globe_crs: CRS, utcnow: datetime) -> np.ndarray:
        raise NotImplementedError()


@dataclass(frozen=True, slots=True)
class NoonPrayer(Prayer):
    def isochrone(self, globe_crs: CRS, utcnow: datetime) -> np.ndarray:
        return isochrone_from_noon_angle(
            globe_crs=globe_crs, utcnow=utcnow, angle=0,
        )


@dataclass(frozen=True, slots=True)
class RefractionPrayer(Prayer):
    angle: float  # degrees after sunrise or before sunset
    pm: bool      # true for any time after noon

    def isochrone(self, globe_crs: CRS, utcnow: datetime) -> np.ndarray:
        """
        :param globe_crs: The coordinate reference system of the globe, used when translating to the
                          night-rotated coordinate system. Typically Geodetic.
        :param utcnow: Timezone-aware datetime used to locate the sun
        :return: A 2*n array of x, y coordinates in degrees; the boundary of "night" where the
                 darkness level is controlled by self.angle.
        """
        night = Nightshade(date=utcnow, delta=2, refraction=self.angle)
        geom, = night.geometries()
        xarray: array
        yarray: array
        xarray, yarray = geom.boundary.coords.xy
        xy = globe_crs.transform_points(
            x=np.frombuffer(xarray),
            y=np.frombuffer(yarray), src_crs=night.crs,
        ).T[:-1]

        # This is wasteful - we always throw away half of the geometry, even if it's symmetrical.
        # However, it's either that or copy a bunch of code from cartopy.
        noon_idx = xy.shape[1]//2
        if self.pm:
            xy = xy[:, :noon_idx]
        else:
            xy = xy[:, noon_idx:]
        return xy


@dataclass(frozen=True, slots=True)
class ShadowPrayer(Prayer):
    shadow: float

    def isochrone(self, globe_crs: CRS, utcnow: datetime) -> np.ndarray:
        sun = SolarPosition.from_time(utcnow=utcnow)
        sun.test()
        A = sun.shadow_angle(shadow=self.shadow)
        return isochrone_from_noon_angle(
            globe_crs=globe_crs, utcnow=utcnow, angle=np.rad2deg(A),
        )


# These definitions can vary significantly; see e.g.
# http://www.praytimes.org/calculation#Fajr_and_Isha
# https://radhifadlillah.com/blog/2020-09-06-calculating-prayer-times/
# 15/15 used in North America by ISNA.
PRAYERS = (
    RefractionPrayer(name='Fajr', colour='orange', angle=-15, pm=False),
    NoonPrayer(name='Dhuhr', colour='yellow'),
    ShadowPrayer(name='Asr', colour='fuchsia', shadow=1),
    RefractionPrayer(name='Maghrib', colour='purple', angle=-0.8333, pm=True),
    RefractionPrayer(name='Isha', colour='blue', angle=-15, pm=True),
)
