import ast
import importlib.util
import inspect
import json
from array import array
from dataclasses import dataclass
from datetime import datetime, timezone, tzinfo
from typing import NamedTuple, Protocol
from zoneinfo import ZoneInfo

import numpy as np
from cartopy.crs import CRS, Geodetic, Orthographic
from cartopy.feature.nightshade import Nightshade, _julian_day, _solar_position
from cartopy.geodesic import Geodesic
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib import pyplot as plt
from matplotlib.patches import Arc
from shapely import Point, Polygon

# Location and time zone of the Kaaba in Mecca
KAABA_COORD = (39.826167, 21.4225)
KAABA_TIMEZONE = ZoneInfo('Asia/Riyadh')

# Day/night graph colours
GEODETIC_COLOUR: tuple[str, str] = ('green', 'lightgreen')
FEATURE_COLOUR: tuple[str, str] = ('black', 'white')
ANNOTATE_COLOUR: tuple[str, str] = ('black', 'white')
HEADING_COLOUR: tuple[str, str] = ('red', 'red')


def check_contains_night(coord: tuple[float, float], night: Nightshade) -> bool:
    geom: Polygon
    geom, = night.geometries()
    point = Point(coord)
    return geom.contains(point)


def time_colour(is_night: bool, colours: tuple[str, str]) -> str:
    """Based on the time, choose light or dark colours for contrast."""
    day, night = colours
    return night if is_night else day


def load_home() -> tuple[float, float]:
    """
    Load the prayer coordinate
    """
    with open('.home.json') as f:
        coords = json.load(f)
    return coords['lon'], coords['lat']


@dataclass(frozen=True)
class Prayer:
    """One salah prayer definition"""
    name: str     # prayer name in romanized Arabic
    colour: str   # time-invariant graphing colour

    def isochrone(self, globe_crs: CRS, utcnow: datetime) -> np.ndarray:
        raise NotImplementedError()


@dataclass(frozen=True, slots=True)
class AnglePrayer(Prayer):
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


class SolarPositionCall(Protocol):
    def __call__(self, date: datetime) -> dict[str, float]:
        ...


def patch_solar_position() -> SolarPositionCall:
    tree = compile(
        source=inspect.getsource(_solar_position),
        filename='nightshade.patch.py',
        mode='exec', flags=ast.PyCF_ONLY_AST,
        dont_inherit=1, optimize=-1,
    )

    class ReturnAll(ast.NodeTransformer):
        def visit_Return(self, node: ast.Return) -> ast.stmt:
            return ast.Return(
                value=ast.Call(
                    func=ast.Name(id='locals', ctx=ast.Load()),
                    args=[], keywords=[],
                ),
            )

    code = compile(
        source=ast.fix_missing_locations(
            ReturnAll().visit(tree)),
        filename='nightshade_patch.py',
        mode='exec', flags=0,
        dont_inherit=1, optimize=-1,
    )
    spec = importlib.util.spec_from_loader(name='nightshade_patch', loader=None)
    nightshade_patch = importlib.util.module_from_spec(spec)
    inner_globals = nightshade_patch.__dict__
    inner_globals.update({
        '__builtins__': {'locals': locals},
        '_julian_day': _julian_day,
        'np': np,
    })
    inner_locals = {}
    exec(code, inner_globals, inner_locals)
    return inner_locals['_solar_position']


solar_position = patch_solar_position()


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
    numerator: float        # for RA
    denominator: float      # for RA
    alpha_sun: float        # RA (right ascension)
    lon: float              # opposite of Greenwich Hour Angle (GHA)

    @classmethod
    def from_time(cls, utcnow: datetime) -> 'SolarPosition':
        return cls(**solar_position(date=utcnow))


@dataclass(frozen=True, slots=True)
class AsrPrayer(Prayer):
    shadow: float

    def isochrone(self, globe_crs: CRS, utcnow: datetime) -> np.ndarray:
        sun = SolarPosition.from_time(utcnow=utcnow)

        # www.praytimes.org/calculation#Asr
        A = np.arccos(
            (
                np.sin(
                    np.arctan(  # arccot(1/x) = arctan(x)
                        1/(
                            self.shadow + np.tan(sun.lambda_ecliptic - sun.delta_sun)
                        )
                    )
                )
                - np.sin(sun.lambda_ecliptic) * np.sin(sun.delta_sun)
            ) / (
                np.cos(sun.lambda_ecliptic) * np.cos(sun.delta_sun)
            )
        )
        return


# These definitions can vary significantly; see e.g.
# http://www.praytimes.org/calculation#Fajr_and_Isha
PRAYERS = (
    AnglePrayer(name='Fajr' ,   colour='orange', angle=-18, pm=False),
    AnglePrayer(name='Dhuhr',   colour='yellow', angle=+90, pm=True),
      AsrPrayer(name='Asr',     colour='lightblue', shadow=1),
    AnglePrayer(name='Maghrib', colour='purple', angle=-0.833, pm=True),
    AnglePrayer(name='Isha',    colour='blue'  , angle=-18, pm=True),
)


class Hemisphere(NamedTuple):
    """
    One hemisphere for graphing. Assumed spherical and not oblate because the orthographic
    projection only accepts spherical. Printed figures like the inverse geodesic still use the more
    accurate elliptical WGS-84.
    """

    name: str                      # Hemisphere name, English or romanized Arabic
    coord: tuple[float, float]     # "here" lon, lat in degrees; hemisphere centred on this
    endpoint: tuple[float, float]  # "there" lon, lat in degrees; geodesic heads there
    crs: Orthographic              # Projective coordinate system for graphing
    geodetic: Geodetic             # Globe coordinate system; this one is spherical
    ax: GeoAxes                    # Plot on this
    timezone: tzinfo | None        # None means local timezone. Used for date display.

    @classmethod
    def make(
        cls,
        name: str,
        index: int,
        coord: tuple[float, float],
        endpoint: tuple[float, float],
        figure: plt.Figure,
        local_timezone: tzinfo | None,
        n_axes: int = 2,
        geodetic: Geodetic | None = None,
    ) -> 'Hemisphere':
        """
        :param index: of the axes within the figure.
        :param figure: to which the new axes will be added
        :param n_axes: in the figure. Typically two, for two hemispheres.
        :param geodetic: if None, then a spherical geodetic will be made.
        :return: A new Hemisphere.
        """
        crs = Orthographic(
            central_longitude=coord[0],
            central_latitude=coord[1], globe=geodetic and geodetic.globe,
        )
        ax = figure.add_subplot(1, n_axes, index, projection=crs)
        return cls(
            name=name, coord=coord, endpoint=endpoint, crs=crs, ax=ax, timezone=local_timezone,
            geodetic=geodetic or Geodetic(globe=crs.globe))

    def plot(
        self,
        dusk: Nightshade,
        night: Nightshade,
        utcnow: datetime,
        include_heading: bool = False,
    ) -> None:
        local_now = utcnow.astimezone(self.timezone)
        is_night = check_contains_night(night=dusk, coord=self.coord)
        self.plot_common(dusk=dusk, night=night, local_now=local_now, is_night=is_night)
        self.plot_prayer_isochrones(utcnow=utcnow)
        if include_heading:
            self.plot_heading(is_night=is_night)

    def plot_common(
        self,
        dusk: Nightshade,
        night: Nightshade,
        local_now: datetime,
        is_night: bool,
    ) -> None:
        """
        Plot the common elements: the surface bitmap, night shading, the geodetic between the prayer
        location and the Kaaba, and the time at the centre.
        :param dusk: to shade between sunset/sunrise and dawn/dusk; low refractive correction
        :param night: to shade full darkness; high refractive correction
        """

        self.ax.set_title(f'{self.name} hemisphere')
        self.ax.stock_img()
        self.ax.add_feature(dusk, zorder=5)
        self.ax.add_feature(night, zorder=5)
        self.ax.gridlines()

        self.ax.plot(
            [self.coord[0], self.endpoint[0]],
            [self.coord[1], self.endpoint[1]],
            transform=self.geodetic, zorder=10, label='Qibla',
            c=time_colour(is_night, GEODETIC_COLOUR),
        )

        self.ax.scatter(
            [self.coord[0]], [self.coord[1]],
            transform=self.geodetic, zorder=11, marker='+',
            c=time_colour(is_night, FEATURE_COLOUR),
        )

        self.ax.text(
            x=self.coord[0], y=self.coord[1], rotation=270,
            s=local_now.strftime('%Y-%m-%d %H:%M:%S %z'),
            transform=self.geodetic, zorder=12, ha='left', va='top',
            color=time_colour(is_night, ANNOTATE_COLOUR),
        )

    def plot_prayer_isochrones(self, utcnow: datetime) -> None:
        """
        Plot all of the prayer isochrone curves. Isochrones are not strictly meridians, don't always
        intersect, and are always partially occluded by the globe. They don't vary in colour based
        on time.
        """
        for prayer in PRAYERS:
            xy = prayer.isochrone(globe_crs=self.geodetic, utcnow=utcnow)
            self.ax.plot(
                *xy, transform=self.geodetic, zorder=10,
                label=prayer.name, c=prayer.colour,
            )

    def inverse_geodesic(self) -> tuple[float, float]:
        """
        Calculate the inverse geodesic parameters between the prayer location and the Kaaba. This
        uses WGS-84 and not the spherical CRS - so we use our own Geodesic instead of self.geodesic.
        :return: Distance in metres, and departing angle in degrees counterclockwise from east.
        """
        geodesic = Geodesic()
        (distance, here_azimuth, other_azimuth), = geodesic.inverse(
            points=self.coord, endpoints=self.endpoint,
        )
        return distance, here_azimuth

    def plot_heading(self, is_night: bool) -> None:
        """
        Plot a 'north' indicator, a heading arc from north to the geodesic, and the departing
        heading in degrees
        """
        distance, here_azimuth = self.inverse_geodesic()

        north = 0
        self.ax.plot(
            [self.coord[0], self.coord[0]],
            [self.coord[1], self.coord[1] + 20],
            transform=self.geodetic, zorder=10,
            c=time_colour(is_night, HEADING_COLOUR),
        )
        self.ax.add_patch(Arc(
            xy=self.coord, width=10, height=10,
            transform=self.geodetic, zorder=10,
            theta1=90 - here_azimuth, theta2=90 - north,
            color=time_colour(is_night, HEADING_COLOUR),
        ))

        self.ax.text(
            x=self.coord[0], y=self.coord[1] + 6, s=f'{here_azimuth:.1f}°',
            transform=self.geodetic, zorder=12,
            color=time_colour(is_night, ANNOTATE_COLOUR),
        )
        self.ax.text(
            x=0.92, y=0.82,
            s=f'Geodesic: {distance*1e-3:,.0f} km',
            transform=self.ax.transAxes, zorder=12,
            color=ANNOTATE_COLOUR[1],  # it's always "night" in space
        )

    def plot_legend(self) -> None:
        """Plot a legend of the geodesic and prayer names"""
        self.ax.legend(
            loc=(1, 0.1), bbox_transform=self.ax.transAxes,
        )


def plot_spherical(
    home_coord: tuple[float, float],
    utcnow: datetime,
) -> plt.Figure:
    """
    Plot two hemispheres in the spherical coordinate system, the left centred on 'home' (the prayer
    location) and the right centred on the Kaaba. These are expected to partially overlap.
    :param home_coord: Home lon, lat in degrees
    :param utcnow: Universal, tz-aware 'now' timestamp
    :return: the plot figure.
    """
    fig: plt.Figure = plt.figure()

    home_hemi = Hemisphere.make(
        name='Home', index=1, coord=home_coord, endpoint=KAABA_COORD,
        figure=fig, local_timezone=None)
    kaaba_hemi = Hemisphere.make(
        name='Kaaba', index=2, coord=KAABA_COORD, endpoint=home_coord,
        figure=fig, local_timezone=KAABA_TIMEZONE, geodetic=home_hemi.geodetic)

    dusk = Nightshade(date=utcnow, delta=2, refraction=0, alpha=0.33)
    night = Nightshade(date=utcnow, delta=2, refraction=-18, alpha=0.33)
    kaaba_hemi.plot(utcnow=utcnow, dusk=dusk, night=night)
    home_hemi.plot(utcnow=utcnow, dusk=dusk, night=night, include_heading=True)
    home_hemi.plot_legend()

    return fig


def main() -> None:
    utcnow = datetime.now().astimezone(timezone.utc)
    home_coord = load_home()

    plt.style.use('dark_background')
    plot_spherical(home_coord=home_coord, utcnow=utcnow)
    plt.show()


if __name__ == '__main__':
    main()
