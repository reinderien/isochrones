import typing
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from .astro import shadow_angle

if typing.TYPE_CHECKING:
    from cartopy.crs import CRS
    from .astro import SolarPosition
    from .types import Degree, DegArray, RadArray


@dataclass(frozen=True)
class Prayer(ABC):
    """One salah prayer definition"""
    name: str     # prayer name in romanized Arabic
    colour: str   # for graphing

    @abstractmethod
    def isochrone(self, globe_crs: 'CRS', sun: 'SolarPosition') -> 'DegArray':
        ...


@dataclass(frozen=True, slots=True)
class NoonPrayer(Prayer):
    @staticmethod
    def make_lon(y: 'RadArray') -> 'RadArray':
        return np.zeros_like(y)

    def isochrone(self, globe_crs: 'CRS', sun: 'SolarPosition') -> 'DegArray':
        return sun.isochrone_from_noon_angle(
            globe_crs=globe_crs, make_lon=self.make_lon,
        )


@dataclass(frozen=True, slots=True)
class RefractionPrayer(Prayer):
    angle: 'Degree'  # degrees after sunrise or before sunset
    pm: bool         # true for any time after noon

    def isochrone(self, globe_crs: 'CRS', sun: 'SolarPosition') -> 'DegArray':
        """
        :param globe_crs: The coordinate reference system of the globe, used when translating to the
                          night-rotated coordinate system. Typically Geodetic.
        :return: A 2*n array of x, y coordinates in degrees; the boundary of "night" where the
                 darkness level is controlled by self.angle.
        """
        return sun.isochrone_from_refraction(
            globe_crs=globe_crs, angle=self.angle, pm=self.pm,
        )


@dataclass(frozen=True, slots=True)
class ShadowPrayer(Prayer):
    shadow: float

    def make_lon(self, y: 'RadArray') -> 'RadArray':
        return shadow_angle(shadow=self.shadow, y=y)

    def isochrone(self, globe_crs: 'CRS', sun: 'SolarPosition') -> 'DegArray':
        return sun.isochrone_from_noon_angle(
            globe_crs=globe_crs, make_lon=self.make_lon,
        )


# These definitions can vary significantly; see e.g.
# http://www.praytimes.org/calculation#Fajr_and_Isha
# https://radhifadlillah.com/blog/2020-09-06-calculating-prayer-times/
# 15/15 used in North America by ISNA.
PRAYERS = (
    RefractionPrayer(name='Fajr', colour='orange', angle=-15, pm=False),
    NoonPrayer(name='Dhuhr', colour='yellow'),
    # "Shafi'i" variant
    ShadowPrayer(name='Asr', colour='fuchsia', shadow=1),
    RefractionPrayer(name='Maghrib', colour='purple', angle=-0.8333, pm=True),
    RefractionPrayer(name='Isha', colour='blue', angle=-15, pm=True),
)
