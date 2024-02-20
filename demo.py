import json
from array import array
from datetime import datetime, timezone, tzinfo
from typing import NamedTuple
from zoneinfo import ZoneInfo

import numpy as np
from cartopy.crs import CRS, Geodetic, Orthographic
from cartopy.feature.nightshade import Nightshade
from cartopy.geodesic import Geodesic
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib import pyplot as plt
from matplotlib.patches import Arc

# Location and time zone of the Kaaba in Mecca
KAABA_COORD = (39.826167, 21.4225)
KAABA_TIMEZONE = ZoneInfo('Asia/Riyadh')

# Night/day graph colours
GEODETIC_COLOUR: tuple[str, str] = ('green', 'lightgreen')
FEATURE_COLOUR: tuple[str, str] = ('black', 'white')
ANNOTATE_COLOUR: tuple[str, str] = ('black', 'white')
HEADING_COLOUR: tuple[str, str] = ('red', 'red')


def time_colour(local_now: datetime, colours: tuple[str, str]) -> str:
    """
    Based on the time, choose light or dark colours for contrast.
    Hacky and hard-coded; could actually pay attention to night data instead.
    """
    day, night = colours
    if 6 <= local_now.hour < 18:
        return day
    return night


def load_home() -> tuple[float, float]:
    """
    Load the prayer coordinate
    """
    with open('.home.json') as f:
        coords = json.load(f)
    return coords['lon'], coords['lat']


class Prayer(NamedTuple):
    """One salah prayer definition"""

    name: str     # prayer name in romanized Arabic
    angle: float  # degrees after sunrise or before sunset
    pm: bool      # true for any time after noon
    colour: str   # time-invariant graphing colour

    def isochrone(self, globe_crs: CRS, utcnow: datetime) -> np.ndarray:
        """
        :param globe_crs: The coordinate reference system of the globe, used when translating to the
                          night-rotated coordinate system. Typically Geodetic.
        :param utcnow: Timezone-aware datetime used to locate the sun
        :return: A 2*n array of x, y coordinates in degrees; the boundary of "night" where the
                 darkness level is controlled by self.angle.
        """
        night = Nightshade(refraction=self.angle, date=utcnow, delta=2)
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


# These definitions can vary significantly; see e.g.
# http://www.praytimes.org/calculation#Fajr_and_Isha
PRAYERS = (
    Prayer(name='Fajr', angle=-18, pm=False, colour='orange'),
    Prayer(name='Dhuhr', angle=+90, pm=True, colour='yellow'),
    Prayer(name='Maghrib', angle=-0.833, pm=True, colour='purple'),
    Prayer(name='Isha', angle=-18, pm=True, colour='blue'),
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

    def plot_common(
        self,
        dusk: Nightshade,
        night: Nightshade,
        utcnow: datetime,
    ) -> None:
        """
        Plot the common elements: the surface bitmap, night shading, the geodetic between the prayer
        location and the kaaba, and the time at the centre.
        :param dusk: to shade between sunset/sunrise and dawn/dusk; low refractive correction
        :param night: to shade full darkness; high refractive correction
        """
        local_now = utcnow.astimezone(self.timezone)

        self.ax.set_title(f'{self.name} hemisphere')
        self.ax.stock_img()
        self.ax.add_feature(dusk, zorder=5)
        self.ax.add_feature(night, zorder=5)
        self.ax.gridlines()

        self.ax.plot(
            [self.coord[0], self.endpoint[0]],
            [self.coord[1], self.endpoint[1]],
            transform=self.geodetic, zorder=10, label='Qibla',
            c=self.time_colour(local_now, GEODETIC_COLOUR),
        )

        self.ax.scatter(
            [self.coord[0]], [self.coord[1]],
            transform=self.geodetic, zorder=11, marker='+',
            c=self.time_colour(local_now, FEATURE_COLOUR),
        )

        self.ax.text(
            x=self.coord[0], y=self.coord[1], rotation=270,
            s=local_now.strftime('%Y-%m-%d %H:%M:%S %z'),
            transform=self.geodetic, zorder=12, ha='left', va='top',
            color=self.time_colour(local_now, ANNOTATE_COLOUR),
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

    def plot_heading(self, utcnow: datetime) -> None:
        """
        Plot a 'north' indicator, a heading arc from north to the geodesic, and the departing
        heading in degrees
        """
        local_now = utcnow.astimezone(self.timezone)
        distance, here_azimuth = self.inverse_geodesic()

        north = 0
        self.ax.plot(
            [self.coord[0], self.coord[0]],
            [self.coord[1], self.coord[1] + 20],
            transform=self.geodetic, zorder=10,
            c=self.time_colour(local_now, HEADING_COLOUR),
        )
        self.ax.add_patch(Arc(
            xy=self.coord, width=10, height=10,
            transform=self.geodetic, zorder=10,
            theta1=90 - here_azimuth, theta2=90 - north,
            color=self.time_colour(local_now, HEADING_COLOUR),
        ))

        self.ax.text(
            x=self.coord[0], y=self.coord[1] + 6, s=f'{here_azimuth:.1f}°',
            transform=self.geodetic, zorder=12,
            color=self.time_colour(local_now, ANNOTATE_COLOUR),
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
    kaaba_hemi.plot_common(utcnow=utcnow, dusk=dusk, night=night)
    home_hemi.plot_common(utcnow=utcnow, dusk=dusk, night=night)
    kaaba_hemi.plot_prayer_isochrones(utcnow=utcnow)
    home_hemi.plot_prayer_isochrones(utcnow=utcnow)
    home_hemi.plot_heading(utcnow=utcnow)
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
