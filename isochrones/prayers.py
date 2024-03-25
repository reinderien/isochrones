import typing
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from .astro import SolarPosition, shadow_angle

if typing.TYPE_CHECKING:
    from cartopy.crs import CRS
    from .types import FloatArray


@dataclass(frozen=True)
class Prayer(ABC):
    """One salah prayer definition"""
    name: str     # prayer name in romanized Arabic
    colour: str   # for graphing

    @abstractmethod
    def isochrone(self, globe_crs: 'CRS', sun: SolarPosition) -> 'FloatArray':
        ...


@dataclass(frozen=True, slots=True)
class NoonPrayer(Prayer):
    @staticmethod
    def make_lon(y: 'FloatArray') -> 'FloatArray':
        return np.zeros_like(y)

    def isochrone(self, globe_crs: 'CRS', sun: SolarPosition) -> 'FloatArray':
        return sun.isochrone_from_noon_angle(
            globe_crs=globe_crs, make_lon=self.make_lon,
        )


@dataclass(frozen=True, slots=True)
class RefractionPrayer(Prayer):
    angle: float  # degrees after sunrise or before sunset
    pm: bool      # true for any time after noon

    def isochrone(self, globe_crs: 'CRS', sun: SolarPosition) -> 'FloatArray':
        """
        :param globe_crs: The coordinate reference system of the globe, used when translating to the
                          night-rotated coordinate system. Typically Geodetic.
        :return: A 2*n array of x, y coordinates in degrees; the boundary of "night" where the
                 darkness level is controlled by self.angle.
        """

        # Adaptation of Nightshade(), but with a coordinate system fixup, and only populating one
        # side of daytime
        refraction = np.deg2rad(self.angle)
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
        if not self.pm:
            x = -x
        
        xy: FloatArray = globe_crs.transform_points(
            x=np.rad2deg(x),
            y=np.rad2deg(y),
            src_crs=sun.rotated_pole,
        ).T[:-1]
        return xy


@dataclass(frozen=True, slots=True)
class ShadowPrayer(Prayer):
    shadow: float

    def make_lon(self, y: 'FloatArray') -> 'FloatArray':
        return shadow_angle(shadow=self.shadow, y=y)

    def isochrone(self, globe_crs: 'CRS', sun: SolarPosition) -> 'FloatArray':
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
    ShadowPrayer(name='Asr', colour='fuchsia', shadow=1),
    RefractionPrayer(name='Maghrib', colour='purple', angle=-0.8333, pm=True),
    RefractionPrayer(name='Isha', colour='blue', angle=-15, pm=True),
)
