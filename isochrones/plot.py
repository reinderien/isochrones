from dataclasses import dataclass
from datetime import datetime, tzinfo
from typing import NamedTuple

from cartopy.crs import Geodetic, Orthographic
from cartopy.feature.nightshade import Nightshade
from cartopy.mpl.feature_artist import FeatureArtist
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib import pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Arc
from matplotlib.text import Text

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
    parallel_on_self: bool
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
            parallel_on_self=parallel_on_self, include_heading=include_heading,
        )


@dataclass(slots=True)
class HemispherePlots(NamedTuple):
    data: HemisphereData
    ax: GeoAxes                    # Plot on this
    dusk_art: FeatureArtist
    night_art: FeatureArtist
    qibla_art: Line2D
    origin_art: PathCollection
    origin_time_art: Text
    is_night: bool

    @classmethod
    def make(
        cls,
        data: HemisphereData,
        initial_dusk: Nightshade,
        initial_night: Nightshade,
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

        cls.plot_invariants(data, ax)

        dusk_art, night_art, qibla_art, origin_art, origin_time_art = cls.setup_common(
            initial_dusk=initial_dusk, initial_night=initial_night,
            initial_is_night=,
        )

        return cls(data=data, ax=ax)

    @classmethod
    def plot_invariants(cls, data: HemisphereData, ax: GeoAxes) -> None:
        ax.set_title(f'{data.name} hemisphere')
        ax.stock_img()
        ax.gridlines()

    @classmethod
    def setup_common(
        cls,
        ax: GeoAxes,
        data: HemisphereData,
        initial_dusk: Nightshade,
        initial_night: Nightshade,
        initial_local: datetime,
        initial_is_night: bool,
    ) -> tuple[
        FeatureArtist, FeatureArtist, Line2D, PathCollection, Text,
    ]:
        dusk_art: FeatureArtist = ax.add_feature(initial_dusk, zorder=5)
        night_art: FeatureArtist = ax.add_feature(initial_night, zorder=5)

        qibla_art: Line2D
        qibla_art, = ax.plot(
            [data.coord[0], data.endpoint[0]],
            [data.coord[1], data.endpoint[1]],
            transform=data.geodetic, zorder=10, label='Qibla',
            c=time_colour(initial_is_night, GEODETIC_COLOUR),
        )

        origin_art: PathCollection = ax.scatter(
            [data.coord[0]], [data.coord[1]],
            transform=data.geodetic, zorder=11, marker='+',
            c=time_colour(initial_is_night, FEATURE_COLOUR),
        )

        origin_time_art: Text = ax.text(
            x=data.coord[0], y=data.coord[1], rotation=270,
            s=initial_local.strftime('%Y-%m-%d %H:%M:%S %z'),
            transform=data.geodetic, zorder=12, ha='left', va='top',
            color=time_colour(initial_is_night, ANNOTATE_COLOUR),
        )

        return dusk_art, night_art, qibla_art, origin_art, origin_time_art

    def update_common(
        self,
        dusk: Nightshade,
        night: Nightshade,
        local_now: datetime,
        is_night: bool,
    ) -> tuple[plt.Artist, ...]:
        self.dusk_art._feature = dusk
        self.night_art._feature = night

        if self.is_night != is_night:
            self.qibla_art.set_color(time_colour(is_night, GEODETIC_COLOUR))
            self.origin_art.set_color(time_colour(is_night, FEATURE_COLOUR))
            self.origin_time_art



        return (
            dusk_art, night_art, qibla_art, origin_art, origin_time_art,
        )

    def update_ecliptic_parallel(self, utcnow: datetime) -> tuple[plt.Artist, ...]:
        xy = ecliptic_parallel(globe_crs=self.data.geodetic, utcnow=utcnow, home=self.data.home)
        parallel_art: Line2D
        parallel_art, = self.ax.plot(
            *xy, transform=self.data.geodetic, zorder=9,
            # label='ecliptic parallel',  # long and not particularly necessary
            c='grey',
        )
        return parallel_art,

    def update_prayer_isochrones(self, utcnow: datetime) -> tuple[plt.Artist, ...]:
        """
        Plot all of the prayer isochrone curves. Isochrones are not strictly meridians, don't always
        intersect, and are always partially occluded by the globe. They don't vary in colour based
        on time.
        """
        artists = []
        for prayer in PRAYERS:
            xy = prayer.isochrone(globe_crs=self.data.geodetic, utcnow=utcnow)
            prayer_art, = self.ax.plot(
                *xy, transform=self.data.geodetic, zorder=10,
                label=prayer.name, c=prayer.colour,
            )
            artists.append(prayer_art)
        return tuple(artists)

    def setup_heading(self) -> None:


    def update_heading(self, is_night: bool) -> tuple[plt.Artist, ...]:
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
