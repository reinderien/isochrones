from array import array
from dataclasses import dataclass
from datetime import datetime

import numpy as np
from cartopy.crs import CRS
from cartopy.feature.nightshade import Nightshade

from .astro import SolarPosition


@dataclass(frozen=True)
class Prayer:
    """One salah prayer definition"""
    name: str     # prayer name in romanized Arabic
    colour: str   # time-invariant graphing colour

    def isochrone(self, globe_crs: CRS, utcnow: datetime) -> np.ndarray:
        raise NotImplementedError()


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
        y = np.linspace(-90, 90, 181)
        xyz = globe_crs.transform_points(
            x=np.full_like(a=y, fill_value=180 + np.rad2deg(A)),
            y=y, src_crs=sun.rotated_pole,
        ).T
        return xyz[:-1]


# These definitions can vary significantly; see e.g.
# http://www.praytimes.org/calculation#Fajr_and_Isha
PRAYERS = (
    RefractionPrayer(name='Fajr', colour='orange', angle=-18, pm=False),
    RefractionPrayer(name='Dhuhr', colour='yellow', angle=+90, pm=True),
    ShadowPrayer(name='Asr', colour='fuchsia', shadow=1),
    RefractionPrayer(name='Maghrib', colour='purple', angle=-0.833, pm=True),
    RefractionPrayer(name='Isha', colour='blue', angle=-18, pm=True),
)
