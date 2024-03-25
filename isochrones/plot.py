import itertools
import logging
import typing
from datetime import datetime, timedelta, timezone, tzinfo

from cartopy.crs import Geodetic, Orthographic
from cartopy.feature import ShapelyFeature
from cartopy.feature.nightshade import Nightshade
from cartopy.mpl.feature_artist import FeatureArtist
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.ticker import MultipleLocator
from matplotlib.transforms import Affine2D

from .astro import (
    KAABA_COORD, KAABA_TIMEZONE, inverse_geodesic, inverse_geodesic_hack, SolarPosition,
)
from .prayers import PRAYERS

if typing.TYPE_CHECKING:
    from cartopy.mpl.geoaxes import GeoAxes
    from matplotlib.collections import PathCollection
    from matplotlib.typing import ColourType
    from pyproj import CRS
    from .types import Coord, FloatArray


GEODESIC_COLOUR = 'green'
FEATURE_COLOUR = 'white'
ANNOTATE_COLOUR = 'white'
HEADING_COLOUR = 'red'

logger = logging.getLogger(__name__)


def hires_arc(
    centre: 'Coord',
    radius: float,
    theta1: float,
    theta2: float,
    n: int | None = None,
    colour: typing.Optional['ColourType'] = None,
    transform: typing.Optional['CRS'] = None,
    zorder: int | None = None,
) -> PathPatch:
    """
    matplotlib's Arc patch looks pretty bad; it has a low and non-configurable resolution. This does
    basically the same thing but with controllable resolution.

    :param centre: x, y pair, centre of circle
    :param radius: from centre
    :param theta1: start angle, degrees counterclockwise from east
    :param theta2: end anglee, degrees counterclockwise from east
    :param n: number of vertices
    :param colour: passed to the patch
    :param transform: passed to the patch
    :param zorder: passed to the patch
    """
    if n is None:
        n = abs(int(round(0.5*(theta2 - theta1))))

    scale = (
        Affine2D()
        .scale(radius)
        .translate(tx=centre[0], ty=centre[1])
    )
    path = Path.arc(
        theta1=theta1, theta2=theta2, n=n,
    ).transformed(scale)
    return PathPatch(
        path=path, edgecolor=colour, facecolor='none', zorder=zorder, transform=transform,
    )


class HemisphereData(typing.NamedTuple):
    """
    One hemisphere's worth of invariant data for graphing. Most graphical elements assume the use of
    a spherical and not ellipsoid CRS because the orthographic projection only accepts spherical.
    Printed figures like the inverse geodesic still use the more accurate elliptical WGS-84.
    """

    name: str            # Hemisphere name, English or romanized Arabic
    coord: 'Coord'       # "here" lon, lat in degrees; hemisphere centred on this
    endpoint: 'Coord'    # "there" lon, lat in degrees; geodesic heads there
    home: 'Coord'        # "home" lon, lat in degrees; ecliptic parallel here
    crs: Orthographic    # Projective coordinate system for graphing
    sphere: Geodetic     # Spherical globe coordinate system
    ellipsoid: Geodetic  # Ellipsoid globe coordinate system

    # Geodesic from coord to endpoint (forward), and reverse (back).
    # All azimuths are in degrees counterclockwise from east. Sphericals are used for plotting;
    # ellipsoidals are used for text.
    geodesic_azm_ell_fwd: float  # Ellipsoid forward
    geodesic_azm_ell_bck: float  # Ellipsoid backward
    geodesic_azm_sph_fwd: float  # Spherical forward
    geodesic_azm_sph_bck: float  # Spherical backward
    geodesic_distance: float  # from coord to endpoint, m

    timezone: tzinfo | None  # None means local timezone. Used for date display.
    is_home: bool            # Is this the home hemisphere?
    include_heading: bool    # Do we include the geodesic heading on this plot?

    @classmethod
    def make(
        cls,
        name: str,
        coord: 'Coord', endpoint: 'Coord',
        local_timezone: tzinfo | None,
        ellipsoid: Geodetic, sphere: Geodetic | None = None,
        geodesic_azm_ell_fwd: float | None = None, geodesic_azm_sph_fwd: float | None = None,
        geodesic_azm_ell_bck: float | None = None, geodesic_azm_sph_bck: float | None = None,
        geodesic_distance: float | None = None,
        is_home: bool = False, include_heading: bool = False,
    ) -> 'HemisphereData':
        """
        :param sphere: if None, then a spherical geodetic will be made.
        :return: A new HemisphereData, invariant; no artist objects.
        """
        crs = Orthographic(
            central_longitude=coord[0],
            central_latitude=coord[1], globe=sphere and sphere.globe,
        )
        sphere = sphere or Geodetic(globe=crs.globe)

        if (
            geodesic_distance is None
            or geodesic_azm_ell_fwd is None or geodesic_azm_sph_fwd is None
            or geodesic_azm_ell_bck is None or geodesic_azm_sph_bck is None
        ):
            # More accurate ellipsoid quantities, for text
            geodesic_azm_ell_fwd, geodesic_azm_ell_bck, geodesic_distance = inverse_geodesic(
                geodetic=ellipsoid, coord=coord, endpoint=endpoint,
            )
            # Less accurate spherical quantities, for display on the orthographic
            # geodesic_azm_sph_fwd, geodesic_azm_sph_bck, _ = inverse_geodesic(
            #     geodetic=sphere, coord=coord, endpoint=endpoint,
            # )
            geodesic_azm_sph_fwd, geodesic_azm_sph_bck = inverse_geodesic_hack(
                sphere=sphere, coord=coord, endpoint=endpoint,
            )

        return cls(
            name=name, coord=coord, endpoint=endpoint, crs=crs, timezone=local_timezone,
            sphere=sphere, ellipsoid=ellipsoid,
            geodesic_azm_ell_fwd=geodesic_azm_ell_fwd, geodesic_azm_ell_bck=geodesic_azm_ell_bck,
            geodesic_azm_sph_fwd=geodesic_azm_sph_fwd, geodesic_azm_sph_bck=geodesic_azm_sph_bck,
            geodesic_distance=geodesic_distance,
            home=coord if is_home else endpoint,
            is_home=is_home, include_heading=include_heading,
        )


class FrameData(typing.NamedTuple):
    utcnow: datetime    # Time used to compute isochrones
    sun: SolarPosition  # Solar coordinates used to compute isochrones
    dusk: Nightshade    # Outer, lighter shading at sunrise/sunset; no refractive correction
    night: Nightshade   # Inner, darker shading at dawn/dusk; high refractive correction
    prayers: tuple['FloatArray', ...]

    @classmethod
    def make(
        cls, utcnow: datetime, sphere: Geodetic,
    ) -> 'FrameData':
        sun = SolarPosition.from_time(utcnow)
        sun.test()

        # this angle should be made to match the angles in the prayer database
        night_angle = 15

        return cls(
            utcnow=utcnow,
            sun=sun,
            dusk=Nightshade(date=utcnow, delta=2, refraction=0, alpha=0.33),
            night=Nightshade(date=utcnow, delta=2, refraction=-night_angle, alpha=0.33),
            prayers=tuple(
                prayer.isochrone(globe_crs=sphere, sun=sun)
                for prayer in PRAYERS
            ),
        )


class HemispherePlots(typing.NamedTuple):
    data: HemisphereData  # Hemisphere invariants

    # Axes and artists that never get reassigned between frames, only updated.
    ax: 'GeoAxes'                 # Plot on this
    dusk_art: FeatureArtist       # Dusk (outer, lighter) shading
    night_art: FeatureArtist      # Night (inner, darker) shading
    geodesic_art: plt.Line2D      # Geodesic "qibla" line, straight due to orthographic projn
    origin_art: 'PathCollection'  # Cross at "here"
    origin_time_art: plt.Text     # Time text at "here"
    prayer_art: tuple[plt.Line2D, ...]  # All prayer isochrones
    north_arc_art: PathPatch | None  # Arc from north to geodesic
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
        :param data: hemisphere invariant geodata
        :param figure: to which the new axes will be added
        :param index: of the axes within the figure.
        :param n_axes: in the figure. Typically two, for two hemispheres.
        :return: A HemispherePlots ready for either plotting or animation.
        """
        logger.info('Setting up artists (%s)...', data.name)

        nrows = 1
        ncols = n_axes
        ax = figure.add_subplot(nrows, ncols, index, projection=data.crs)

        # This dominates the load time
        cls.plot_invariants(data, ax)

        (
            dusk_art, night_art, geodesic_art, origin_art, origin_time_art,
        ) = cls.setup_common(ax=ax, data=data)
        prayer_art = cls.setup_prayer_isochrones(ax=ax, data=data)
        north_arc_art, heading_art = (
            cls.setup_heading(ax=ax, data=data)
            if data.include_heading else (None, None)
        )

        return cls(
            data=data, ax=ax,
            dusk_art=dusk_art, night_art=night_art, geodesic_art=geodesic_art,
            origin_art=origin_art, origin_time_art=origin_time_art,
            prayer_art=prayer_art, north_arc_art=north_arc_art,
            heading_art=heading_art,
        )

    def update(self, data: FrameData) -> tuple[plt.Artist, ...]:
        """
        Update all artists as necessary, for either a static plot or one frame of an animation.
        :return: All changed artists.
        """
        changed = (
            self.update_common(data)
            + self.update_prayer_isochrones(data)
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
        ax.gridlines(
            xlocs=MultipleLocator(base=45),
            ylocs=MultipleLocator(base=30), color='white', alpha=0.4,
        )
        ax.gridlines(
            xlocs=data.home[:1],
            ylocs=data.home[1:], color='white', alpha=0.8,
        )

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
            transform=data.sphere, zorder=10, label='Qibla',
            c=GEODESIC_COLOUR,
        )

        origin_art: plt.PathCollection = ax.scatter(
            [data.coord[0]], [data.coord[1]],
            transform=data.sphere, zorder=11, marker='+',
            c=FEATURE_COLOUR,
        )

        origin_time_art: plt.Text = ax.text(
            x=data.coord[0], y=data.coord[1], rotation=-60, s='',
            transform=data.sphere, zorder=12, ha='left', va='top',
            c=ANNOTATE_COLOUR,
        )

        return dusk_art, night_art, geodesic_art, origin_art, origin_time_art

    def update_common(self, data: FrameData) -> tuple[plt.Artist, ...]:
        self.dusk_art._feature = data.dusk
        self.night_art._feature = data.night

        local_now = data.utcnow.astimezone(self.data.timezone)
        self.origin_time_art.set_text(
            local_now.strftime('%Y-%m-%d %H:%M:%S %z')
        )

        # This includes artists that haven't actually updated, but should be redrawn on top of the
        # nightshades
        return (
            self.dusk_art, self.night_art, self.geodesic_art, self.origin_art, self.origin_time_art,
        )

    @classmethod
    def setup_prayer_isochrones(
        cls, ax: 'GeoAxes', data: HemisphereData,
    ) -> tuple[plt.Line2D, ...]:
        """
        Plot all of the prayer isochrone curves. Isochrones are not strictly meridians, don't always
        intersect, and are always partially occluded by the globe. They move over time.
        """
        artists = tuple(itertools.chain.from_iterable(
            ax.plot(
                [], transform=data.sphere, zorder=10,
                label=prayer.name, c=prayer.colour,
            )
            for prayer in PRAYERS
        ))
        return artists

    def update_prayer_isochrones(self, data: FrameData) -> tuple[plt.Artist, ...]:
        for xy, artist in zip(data.prayers, self.prayer_art):
            artist.set_data(*xy)
        return self.prayer_art

    @classmethod
    def setup_heading(
        cls, ax: 'GeoAxes', data: HemisphereData,
    ) -> tuple[
        PathPatch, plt.Text,
    ]:
        """
        Plot a heading arc from north to the geodesic, and the departing
        heading in degrees. All frame-invariant.
        """

        north = 90
        north_arc_art = hires_arc(
            centre=data.coord, radius=10,
            transform=data.sphere, zorder=10,
            theta1=data.geodesic_azm_sph_fwd,
            theta2=north,
            colour=HEADING_COLOUR,
        )
        ax.add_patch(north_arc_art)

        heading_art: plt.Text = ax.text(
            x=data.coord[0], y=data.coord[1] + 6, s=f'{data.geodesic_azm_ell_fwd:.1f}°',
            transform=data.sphere, zorder=12,
            color=ANNOTATE_COLOUR,
        )
        ax.text(
            x=0.92, y=0.82,
            s=f'Geodesic: {data.geodesic_distance*1e-3:,.0f} km',
            transform=ax.transAxes, zorder=12,
            color=ANNOTATE_COLOUR,  # it's always "night" in space
        )
        return north_arc_art, heading_art

    def update_heading(self) -> tuple[plt.Artist, ...]:
        # None of these have updated, but they need to be redrawn on top of the nightshades
        return tuple(
            artist
            for artist in (self.north_arc_art, self.heading_art)
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
    :return: the plot figure and both plot data objects
    """
    plt.style.use('dark_background')
    fig: plt.Figure = plt.figure()

    ellipsoid = Geodetic()
    home_data = HemisphereData.make(
        name='Home', coord=home_coord, endpoint=KAABA_COORD, local_timezone=None,
        ellipsoid=ellipsoid, include_heading=True, is_home=True,
    )
    kaaba_data = HemisphereData.make(
        name='Kaaba', coord=KAABA_COORD, endpoint=home_coord, local_timezone=KAABA_TIMEZONE,
        ellipsoid=ellipsoid, sphere=home_data.sphere,
        geodesic_azm_ell_fwd=home_data.geodesic_azm_ell_bck,
        geodesic_azm_ell_bck=home_data.geodesic_azm_ell_fwd,
        geodesic_azm_sph_fwd=home_data.geodesic_azm_sph_bck,
        geodesic_azm_sph_bck=home_data.geodesic_azm_sph_fwd,
        geodesic_distance=home_data.geodesic_distance,
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
    data = FrameData.make(utcnow=utcnow, sphere=home_plot.data.sphere)
    return home_plot.update(data) + kaaba_plot.update(data)


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
        track_system = time_factor == 1
    else:
        track_system = False

    def update(frame: int) -> tuple[plt.Artist, ...]:
        if track_system:
            # Don't trust that frame * interval == elapsed... because it isn't
            virtual_time = datetime.now().astimezone(timezone.utc)
        else:
            virtual_time = start_utc + timedelta(
                seconds=time_factor*frame,
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
