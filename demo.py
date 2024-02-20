import json
from array import array
from datetime import datetime, timezone, tzinfo
from typing import NamedTuple
from zoneinfo import ZoneInfo

import numpy as np
from cartopy.crs import Geodetic, Orthographic
from cartopy.feature.nightshade import Nightshade
from cartopy.geodesic import Geodesic
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib import pyplot as plt
from matplotlib.patches import Arc

KAABA_COORD = (39.826167, 21.4225)
KAABA_TIMEZONE = ZoneInfo('Asia/Riyadh')

GEODETIC_COLOUR: tuple[str, str] = ('green', 'lightgreen')
FEATURE_COLOUR: tuple[str, str] = ('black', 'white')
ANNOTATE_COLOUR: tuple[str, str] = ('black', 'white')
HEADING_COLOUR: tuple[str, str] = ('red', 'red')


def load_home() -> tuple[float, float]:
    with open('.home.json') as f:
        coords = json.load(f)
    return coords['lon'], coords['lat']


class Hemisphere(NamedTuple):
    name: str
    coord: tuple[float, float]
    endpoint: tuple[float, float]
    crs: Orthographic
    geodetic: Geodetic
    ax: GeoAxes
    timezone: tzinfo | None

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
        crs = Orthographic(
            central_longitude=coord[0],
            central_latitude=coord[1], globe=geodetic and geodetic.globe,
        )
        ax = figure.add_subplot(1, n_axes, index, projection=crs)
        return cls(
            name=name, coord=coord, endpoint=endpoint, crs=crs, ax=ax, timezone=local_timezone,
            geodetic=geodetic or Geodetic(globe=crs.globe))

    @staticmethod
    def time_colour(local_now: datetime, colours: tuple[str, str]) -> str:
        day, night = colours
        if 6 <= local_now.hour < 18:  # hacky; could actually pay attention to night data
            return day
        return night

    def plot_common(
        self,
        night: Nightshade,
        utcnow: datetime,
    ) -> None:
        local_now = utcnow.astimezone(self.timezone)

        self.ax.set_title(f'{self.name} hemisphere')
        self.ax.stock_img()
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

    def plot_salah_meridians(self, night: Nightshade) -> None:
        geom, = night.geometries()
        xarray: array
        yarray: array
        xarray, yarray = geom.boundary.coords.xy
        xy = self.geodetic.transform_points(
            x=np.frombuffer(xarray),
            y=np.frombuffer(yarray), src_crs=night.crs,
        ).T[:-1]
        noon = xy.shape[1]//2
        self.ax.plot(
            *xy[:, noon:], transform=self.geodetic, zorder=10,
            label='Fajr', c='orange',
        )
        self.ax.plot(
            *self.noon_meridian(night), transform=self.geodetic, zorder=10,
            label='Zuhr', c='yellow',
        )
        self.ax.plot(
            *xy[:, :noon], transform=self.geodetic, zorder=10,
            label='Maghrib', c='purple',
        )

    def inverse_geodesic(self) -> tuple[float, float]:
        # WGS-84, not spherical - so we don't pass in class params
        geodesic = Geodesic()
        (distance, here_azimuth, other_azimuth), = geodesic.inverse(
            points=self.coord, endpoints=self.endpoint,
        )
        return distance, here_azimuth

    def noon_meridian(self, night: Nightshade) -> np.ndarray:
        y = np.linspace(start=-90, stop=90, num=51)
        x = np.full_like(a=y, fill_value=180)
        return self.geodetic.transform_points(
            x=x, y=y, src_crs=night.crs,
        ).T[:-1]

    def plot_heading(self, utcnow: datetime) -> None:
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
        self.ax.legend(
            loc=(1, 0.1), bbox_transform=self.ax.transAxes,
        )


def plot_spherical(
    home_coord: tuple[float, float],
    utcnow: datetime,
) -> plt.Figure:
    fig: plt.Figure = plt.figure()

    home_hemi = Hemisphere.make(
        name='Home', index=1, coord=home_coord, endpoint=KAABA_COORD,
        figure=fig, local_timezone=None)
    kaaba_hemi = Hemisphere.make(
        name='Kaaba', index=2, coord=KAABA_COORD, endpoint=home_coord,
        figure=fig, local_timezone=KAABA_TIMEZONE, geodetic=home_hemi.geodetic)

    night = Nightshade(date=utcnow)
    kaaba_hemi.plot_common(night=night, utcnow=utcnow)
    home_hemi.plot_common(night=night, utcnow=utcnow)
    kaaba_hemi.plot_salah_meridians(night=night)
    home_hemi.plot_salah_meridians(night=night)
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
