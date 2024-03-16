from datetime import datetime, tzinfo
from typing import NamedTuple

from cartopy.crs import Geodetic, Orthographic
from cartopy.feature.nightshade import Nightshade
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib import pyplot as plt
from matplotlib.patches import Arc

from .astro import KAABA_COORD, KAABA_TIMEZONE, check_contains_night, ecliptic_parallel, inverse_geodesic
from .prayers import PRAYERS

# Day/night graph colours
GEODETIC_COLOUR: tuple[str, str] = ('green', 'lightgreen')
FEATURE_COLOUR: tuple[str, str] = ('black', 'white')
ANNOTATE_COLOUR: tuple[str, str] = ('black', 'white')
HEADING_COLOUR: tuple[str, str] = ('red', 'red')


def time_colour(is_night: bool, colours: tuple[str, str]) -> str:
    """Based on the time, choose light or dark colours for contrast."""
    day, night = colours
    return night if is_night else day


class HemisphereData(NamedTuple):
    """
    One hemisphere for graphing. Assumed spherical and not oblate because the orthographic
    projection only accepts spherical. Printed figures like the inverse geodesic still use the more
    accurate elliptical WGS-84.
    """

    name: str                      # Hemisphere name, English or romanized Arabic
    coord: tuple[float, float]     # "here" lon, lat in degrees; hemisphere centred on this
    endpoint: tuple[float, float]  # "there" lon, lat in degrees; geodesic heads there
    home: tuple[float, float]
    crs: Orthographic              # Projective coordinate system for graphing
    geodetic: Geodetic             # Globe coordinate system; this one is spherical
    timezone: tzinfo | None        # None means local timezone. Used for date display.
    include_heading: bool

    @classmethod
    def make(
        cls,
        name: str,
        coord: tuple[float, float],
        endpoint: tuple[float, float],
        local_timezone: tzinfo | None,
        geodetic: Geodetic | None = None,
        parallel_on_self: bool = False,
        include_heading: bool = False,
    ) -> 'HemisphereData':
        """
        :param geodetic: if None, then a spherical geodetic will be made.
        :return: A new Hemisphere.
        """
        crs = Orthographic(
            central_longitude=coord[0],
            central_latitude=coord[1], globe=geodetic and geodetic.globe,
        )
        return cls(
            name=name, coord=coord, endpoint=endpoint, crs=crs, timezone=local_timezone,
            geodetic=geodetic or Geodetic(globe=crs.globe),
            home=coord if parallel_on_self else endpoint,
            include_heading=include_heading,
        )


class HemispherePlots(NamedTuple):
    data: HemisphereData
    ax: GeoAxes                    # Plot on this

    @classmethod
    def make(
        cls,
        data: HemisphereData,
        figure: plt.Figure,
        index: int,
        n_axes: int = 2,
    ) -> 'HemispherePlots':
        """
        :param data:
        :param index: of the axes within the figure.
        :param figure: to which the new axes will be added
        :param n_axes: in the figure. Typically two, for two hemispheres.
        :return:
        """
        ax = figure.add_subplot(nrows=1, ncols=n_axes, index=index, projection=data.crs)

        return cls(
        )

    @classmethod
    def setup_invariants(cls) -> None:
        self.ax.set_title(f'{self.name} hemisphere')
        self.heading()

    #
    #
    # @classmethod
    # def plot(
    #     cls,
    #     # dusk: Nightshade,
    #     # night: Nightshade,
    #     # utcnow: datetime,
    # ) -> None:
    #     # local_now = utcnow.astimezone(self.timezone)
    #     # is_night = check_contains_night(night=dusk, local_now=local_now, coord=self.coord)
    #
    #     self.plot_common(dusk=dusk, night=night, local_now=local_now, is_night=is_night)
    #     self.plot_ecliptic_parallel(home=home, utcnow=utcnow)
    #     self.plot_prayer_isochrones(utcnow=utcnow)
    #     if include_heading:
    #         self.plot_heading(is_night=is_night)
    #
    # def plot_common(
    #     self,
    #     dusk: Nightshade,
    #     night: Nightshade,
    #     local_now: datetime,
    #     is_night: bool,
    # ) -> None:
        """
        Plot the common elements: the surface bitmap, night shading, the geodetic between the prayer
        location and the Kaaba, and the time at the centre.
        :param dusk: to shade between sunset/sunrise and dawn/dusk; low refractive correction
        :param night: to shade full darkness; high refractive correction
        """

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

    def plot_ecliptic_parallel(self, utcnow: datetime, home: tuple[float, float]) -> None:
        xy = ecliptic_parallel(globe_crs=self.geodetic, utcnow=utcnow, home=home)
        self.ax.plot(
            *xy, transform=self.geodetic, zorder=9,
            # label='ecliptic parallel',  # long and not particularly necessary
            c='grey',
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

    def plot_heading(self, is_night: bool) -> None:
        """
        Plot a 'north' indicator, a heading arc from north to the geodesic, and the departing
        heading in degrees
        """
        distance, here_azimuth = inverse_geodesic(coord=self.coord, endpoint=self.endpoint)

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
    plt.style.use('dark_background')
    fig: plt.Figure = plt.figure()

    home_data = HemisphereData.make(
        name='Home', index=1, coord=home_coord, endpoint=KAABA_COORD,
        figure=fig, local_timezone=None)
    kaaba_data = HemisphereData.make(
        name='Kaaba', index=2, coord=KAABA_COORD, endpoint=home_coord,
        figure=fig, local_timezone=KAABA_TIMEZONE, geodetic=home_data.geodetic)

    # this angle should be made to match the angles in the prayer database
    night = Nightshade(date=utcnow, delta=2, refraction=-15, alpha=0.33)
    dusk = Nightshade(date=utcnow, delta=2, refraction=0, alpha=0.33)
    kaaba_hemi.plot(utcnow=utcnow, dusk=dusk, night=night)
    home_hemi.plot(utcnow=utcnow, dusk=dusk, night=night,
                   include_heading=True, parallel_on_self=True)
    home_hemi.plot_legend()

    return fig
