import itertools
import time
import typing
from datetime import datetime, timedelta, timezone, tzinfo
from typing import NamedTuple

from cartopy.crs import Geodetic, Orthographic
from cartopy.feature import ShapelyFeature
from cartopy.feature.nightshade import Nightshade
from cartopy.mpl.feature_artist import FeatureArtist
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Arc

from .astro import (
    KAABA_COORD, KAABA_TIMEZONE, ecliptic_parallel, inverse_geodesic,
)
from .prayers import PRAYERS

if typing.TYPE_CHECKING:
    from cartopy.mpl.geoaxes import GeoAxes
    from matplotlib.collections import PathCollection
    from .types import Coord


GEODESIC_COLOUR = 'green'
FEATURE_COLOUR = 'white'
ANNOTATE_COLOUR = 'white'
HEADING_COLOUR = 'red'


class HemisphereData(NamedTuple):
    """
    One hemisphere's worth of invariant data for graphing. Assumed spherical and not oblate because
    the orthographic projection only accepts spherical. Printed figures like the inverse geodesic
    still use the more accurate elliptical WGS-84.
    """

    name: str                # Hemisphere name, English or romanized Arabic
    coord: 'Coord'           # "here" lon, lat in degrees; hemisphere centred on this
    endpoint: 'Coord'        # "there" lon, lat in degrees; geodesic heads there
    home: 'Coord'            # "home" lon, lat in degrees; ecliptic parallel here
    crs: Orthographic        # Projective coordinate system for graphing
    geodetic: Geodetic       # Globe coordinate system; this one is spherical
    timezone: tzinfo | None  # None means local timezone. Used for date display.
    parallel_on_self: bool   # Is the ecliptic parallel on 'coord'?
    include_heading: bool    # Do we include the geodesic heading on this plot?

    @classmethod
    def make(
        cls,
        name: str,
        coord: 'Coord',
        endpoint: 'Coord',
        local_timezone: tzinfo | None,
        geodetic: Geodetic | None = None,
        parallel_on_self: bool = False,
        include_heading: bool = False,
    ) -> 'HemisphereData':
        """
        :param geodetic: if None, then a spherical geodetic will be made.
        :return: A new HemisphereData, invariant; no artist objects.
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


class HemispherePlots(NamedTuple):
    data: HemisphereData  # Hemisphere invariants

    # Axes and artists that never get reassigned between frames, only updated.
    ax: 'GeoAxes'                 # Plot on this
    dusk_art: FeatureArtist       # Dusk (outer, lighter) shading
    night_art: FeatureArtist      # Night (inner, darker) shading
    geodesic_art: plt.Line2D      # Geodesic "qibla" line, straight due to orthographic projn
    origin_art: 'PathCollection'  # Cross at "here"
    origin_time_art: plt.Text     # Time text at "here"
    parallel_art: plt.Line2D      # Ecliptic parallel through home
    prayer_art: tuple[plt.Line2D, ...]  # All prayer isochrones
    north_line_art: plt.Line2D | None
    north_arc_art: Arc | None     # Arc from north to geodesic
    heading_art: plt.Text | None  # Heading text, if we have it

    @classmethod
    def make(
        cls,
        data: HemisphereData,
        figure: plt.Figure,
        index: int,
        n_axes: int = 2,
    ) -> 'HemispherePlots':
        """
        :param index: of the axes within the figure.
        :param figure: to which the new axes will be added
        :param n_axes: in the figure. Typically two, for two hemispheres.
        :return: A HemispherePlots ready for either plotting or animation.
        """
        nrows = 1
        ncols = n_axes
        ax = figure.add_subplot(nrows, ncols, index, projection=data.crs)

        cls.plot_invariants(data, ax)

        (
            dusk_art, night_art, geodesic_art, origin_art, origin_time_art,
        ) = cls.setup_common(ax=ax, data=data)
        parallel_art = cls.setup_ecliptic_parallel(ax=ax, data=data)
        prayer_art = cls.setup_prayer_isochrones(ax=ax, data=data)
        north_line_art, north_arc_art, heading_art = (
            cls.setup_heading(ax=ax, data=data)
            if data.include_heading else (None, None, None)
        )

        return cls(
            data=data, ax=ax,
            dusk_art=dusk_art, night_art=night_art, geodesic_art=geodesic_art,
            origin_art=origin_art, origin_time_art=origin_time_art, parallel_art=parallel_art,
            prayer_art=prayer_art, north_line_art=north_line_art, north_arc_art=north_arc_art,
            heading_art=heading_art,
        )

    def update(
        self, dusk: Nightshade, night: Nightshade, utcnow: datetime,
    ) -> tuple[plt.Artist, ...]:
        """
        Update all artists as necessary, for either a static plot or one frame of an animation.
        :param dusk: Outer, lighter shading at sunrise/sunset; no refractive correction
        :param night: Inner, darker shading at dawn/dusk; high refractive correction
        :param utcnow: Time used to compute ecliptic parallel and isochrones
        :return: All changed artists.
        """
        local_now = utcnow.astimezone(self.data.timezone)

        changed = (
            self.update_common(local_now=local_now, dusk=dusk, night=night)
            + self.update_ecliptic_parallel(utcnow=utcnow)
            + self.update_prayer_isochrones(utcnow=utcnow)
        )
        if self.data.include_heading:
            changed += self.update_heading()

        return changed

    @classmethod
    def plot_invariants(cls, data: HemisphereData, ax: 'GeoAxes') -> None:
        """
        All plot operations that never change frame-to-frame.
        """
        ax.set_title(f'{data.name} hemisphere')
        ax.background_img(name='blue-marble-next-generation', resolution='high')
        ax.gridlines()

    @classmethod
    def setup_common(
        cls, ax: 'GeoAxes', data: HemisphereData,
    ) -> tuple[
        FeatureArtist, FeatureArtist, plt.Line2D, 'PathCollection', plt.Text,
    ]:
        """
        Plot the common elements: the surface bitmap, night shading, the geodetic between the prayer
        location and the Kaaba, and the time at the centre.
        """
        null_feature = ShapelyFeature(geometries=[], crs=data.crs)
        dusk_art: FeatureArtist = ax.add_feature(null_feature, zorder=5)
        night_art: FeatureArtist = ax.add_feature(null_feature, zorder=5)

        geodesic_art: plt.Line2D
        geodesic_art, = ax.plot(
            [data.coord[0], data.endpoint[0]],
            [data.coord[1], data.endpoint[1]],
            transform=data.geodetic, zorder=10, label='Qibla',
            c=GEODESIC_COLOUR,
        )

        origin_art: plt.PathCollection = ax.scatter(
            [data.coord[0]], [data.coord[1]],
            transform=data.geodetic, zorder=11, marker='+',
            c=FEATURE_COLOUR,
        )

        origin_time_art: plt.Text = ax.text(
            x=data.coord[0], y=data.coord[1], rotation=270, s='',
            transform=data.geodetic, zorder=12, ha='left', va='top',
            c=ANNOTATE_COLOUR,
        )

        return dusk_art, night_art, geodesic_art, origin_art, origin_time_art

    def update_common(
        self, dusk: Nightshade, night: Nightshade, local_now: datetime,
    ) -> tuple[plt.Artist, ...]:
        self.dusk_art._feature = dusk
        self.night_art._feature = night
        self.origin_time_art.set_text(
            local_now.strftime('%Y-%m-%d %H:%M:%S %z')
        )

        # This includes artists that haven't actually updated, but should be redrawn on top of the
        # nightshades
        return (
            self.dusk_art, self.night_art, self.geodesic_art, self.origin_art, self.origin_time_art,
        )

    @classmethod
    def setup_ecliptic_parallel(
        cls, ax: 'GeoAxes', data: HemisphereData,
    ) -> plt.Line2D:
        """
        Plot a parallel in the ecliptic (not rotational) frame, through "home". Depicts the progress
        of solar time across the globe for home's latitude.
        """
        parallel_art: plt.Line2D
        parallel_art, = ax.plot(
            [], transform=data.geodetic, zorder=9,
            # label='ecliptic parallel',  # long and not particularly necessary
            c='grey',
        )
        return parallel_art

    def update_ecliptic_parallel(self, utcnow: datetime) -> tuple[plt.Artist, ...]:
        xy = ecliptic_parallel(globe_crs=self.data.geodetic, utcnow=utcnow, home=self.data.home)
        self.parallel_art.set_data(*xy)
        return self.parallel_art,

    @classmethod
    def setup_prayer_isochrones(
        cls, ax: 'GeoAxes', data: HemisphereData,
    ) -> tuple[plt.Line2D, ...]:
        """
        Plot all of the prayer isochrone curves. Isochrones are not strictly meridians, don't always
        intersect, and are always partially occluded by the globe. They vary over time in position
        but not colour.
        """
        artists = tuple(itertools.chain.from_iterable(
            ax.plot(
                [], transform=data.geodetic, zorder=10,
                label=prayer.name, c=prayer.colour,
            )
            for prayer in PRAYERS
        ))
        return artists

    def update_prayer_isochrones(self, utcnow: datetime) -> tuple[plt.Artist, ...]:
        for prayer, artist in zip(PRAYERS, self.prayer_art):
            xy = prayer.isochrone(globe_crs=self.data.geodetic, utcnow=utcnow)
            artist.set_data(*xy)
        return self.prayer_art

    @classmethod
    def setup_heading(
        cls, ax: 'GeoAxes', data: HemisphereData,
    ) -> tuple[
        plt.Line2D, Arc, plt.Text,
    ]:
        """
        Plot a 'north' indicator, a heading arc from north to the geodesic, and the departing
        heading in degrees. Most of this is frame-invariant except the heading text colour.
        """
        distance, here_azimuth = inverse_geodesic(coord=data.coord, endpoint=data.endpoint)

        north = 0
        north_line_art, = ax.plot(
            [data.coord[0], data.coord[0]],
            [data.coord[1], data.coord[1] + 20],
            transform=data.geodetic, zorder=10,
            c=HEADING_COLOUR,
        )
        north_arc_art = Arc(
            xy=data.coord, width=10, height=10,
            transform=data.geodetic, zorder=10,
            theta1=90 - here_azimuth, theta2=90 - north,
            color=HEADING_COLOUR,
        )
        ax.add_patch(north_arc_art)

        heading_art: plt.Text = ax.text(
            x=data.coord[0], y=data.coord[1] + 6, s=f'{here_azimuth:.1f}°',
            transform=data.geodetic, zorder=12,
            color=ANNOTATE_COLOUR,
        )
        ax.text(
            x=0.92, y=0.82,
            s=f'Geodesic: {distance*1e-3:,.0f} km',
            transform=ax.transAxes, zorder=12,
            color=ANNOTATE_COLOUR,  # it's always "night" in space
        )
        return north_line_art, north_arc_art, heading_art

    def update_heading(self) -> tuple[plt.Artist, ...]:
        # None of these have updated, but they need to be redrawn on top of the nightshades
        return tuple(
            artist
            for artist in (self.north_line_art, self.north_arc_art, self.heading_art)
            if artist is not None
        )

    def plot_legend(self) -> None:
        """
        Plot a legend of the geodesic and prayer names. Assume none of these change colour between
        day and night.
        """
        self.ax.legend(
            loc=(1, 0.1), bbox_transform=self.ax.transAxes,
        )


def setup_spherical(home_coord: 'Coord') -> tuple[
    plt.Figure,
    HemispherePlots,
    HemispherePlots,
]:
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
        name='Home', coord=home_coord, endpoint=KAABA_COORD, local_timezone=None,
        include_heading=True, parallel_on_self=True,
    )
    kaaba_data = HemisphereData.make(
        name='Kaaba', coord=KAABA_COORD, endpoint=home_coord, local_timezone=KAABA_TIMEZONE,
        geodetic=home_data.geodetic,
    )

    home_plot = HemispherePlots.make(figure=fig, index=1, data=home_data)
    kaaba_plot = HemispherePlots.make(figure=fig, index=2, data=kaaba_data)
    home_plot.plot_legend()

    return fig, home_plot, kaaba_plot


def update_spherical(
    home_plot: HemispherePlots,
    kaaba_plot: HemispherePlots,
    utcnow: datetime,
) -> tuple[plt.Artist, ...]:
    """
    Update the figure for one frame
    :param utcnow: Universal, tz-aware 'now' timestamp
    """
    # this angle should be made to match the angles in the prayer database
    night = Nightshade(date=utcnow, delta=2, refraction=-15, alpha=0.33)
    dusk = Nightshade(date=utcnow, delta=2, refraction=0, alpha=0.33)
    return (
        home_plot.update(utcnow=utcnow, dusk=dusk, night=night) +
        kaaba_plot.update(utcnow=utcnow, dusk=dusk, night=night)
    )


def plot_spherical(
    home_coord: 'Coord',
    utcnow: datetime | None = None,
) -> plt.Figure:
    """
    Set up the figure and update it once, for a static plot
    :param home_coord: Home lon, lat in degrees
    :param utcnow: Universal, tz-aware 'now' timestamp
    """
    if utcnow is None:
        utcnow = datetime.now().astimezone(timezone.utc)
    fig, home_plot, kaaba_plot = setup_spherical(home_coord)
    update_spherical(home_plot=home_plot, kaaba_plot=kaaba_plot, utcnow=utcnow)
    return fig


def animate_spherical(
    home_coord: 'Coord',
    start_utc: datetime | None = None,
    time_factor: float = 1,
) -> tuple[plt.Figure, FuncAnimation]:
    fig, home_plot, kaaba_plot = setup_spherical(home_coord)

    if start_utc is None:
        start_utc = datetime.now().astimezone(timezone.utc)
    start = time.monotonic()

    def update(frame: int) -> tuple[plt.Artist, ...]:
        # Don't trust that frame * interval == elapsed... because it isn't

        if time_factor == 1:
            virtual_time = datetime.now().astimezone(timezone.utc)
        else:
            virtual_time = start_utc + timedelta(
                seconds=time_factor*(time.monotonic() - start),
            )
        return update_spherical(
            home_plot=home_plot, kaaba_plot=kaaba_plot,
            utcnow=virtual_time,
        )

    anim = FuncAnimation(
        fig=fig, func=update, repeat=False, blit=True, cache_frame_data=False,
        interval=750,
    )

    return fig, anim
