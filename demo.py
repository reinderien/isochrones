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

    GEODETIC_COLOUR = ('green', 'lightgreen')
    FEATURE_COLOUR = ('black', 'white')
    ANNOTATE_COLOUR = ('black', 'white')
    HEADING_COLOUR = ('red', 'red')

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
            transform=self.geodetic, zorder=10,
            c=self.time_colour(local_now, self.GEODETIC_COLOUR),
        )

        self.ax.scatter(
            [self.coord[0]], [self.coord[1]],
            transform=self.geodetic, zorder=11, marker='+',
            c=self.time_colour(local_now, self.FEATURE_COLOUR),
        )

        self.ax.text(
            x=self.coord[0], y=self.coord[1], rotation=270,
            s=local_now.strftime('%Y-%m-%d %H:%M:%S %z'),
            transform=self.geodetic, zorder=12, ha='left', va='top',
            color=self.time_colour(local_now, self.ANNOTATE_COLOUR),
        )

    def plot_salah_isocurves(self, night: Nightshade) -> None:
        geom, = night.geometries()
        xarray: array
        yarray: array
        xarray, yarray = geom.boundary.coords.xy
        x, y, z = self.geodetic.transform_points(
            x=np.frombuffer(xarray),
            y=np.frombuffer(yarray), src_crs=night.crs,
        ).T
        is_sunset = slice(x.size//2)
        is_sunrise = slice(x.size//2, None)
        self.ax.plot(
            x[is_sunrise], y[is_sunrise], transform=self.geodetic, zorder=10, c='orange',
        )
        self.ax.plot(
            x[is_sunset], y[is_sunset], transform=self.geodetic, zorder=10, c='purple',
        )

    def inverse_geodesic(self) -> tuple[float, float]:
        # WGS-84, not spherical - so we don't pass in class params
        geodesic = Geodesic()
        (distance, here_azimuth, other_azimuth), = geodesic.inverse(
            points=self.coord, endpoints=self.endpoint,
        )
        return distance, here_azimuth

    def plot_heading(self, utcnow: datetime) -> None:
        local_now = utcnow.astimezone(self.timezone)
        distance, here_azimuth = self.inverse_geodesic()

        north = 0
        self.ax.plot(
            [self.coord[0], self.coord[0]],
            [self.coord[1], self.coord[1] + 20],
            transform=self.geodetic, zorder=10,
            c=self.time_colour(local_now, self.HEADING_COLOUR),
        )
        self.ax.add_patch(Arc(
            xy=self.coord, width=10, height=10,
            transform=self.geodetic, zorder=10,
            theta1=90 - here_azimuth, theta2=90 - north,
            color=self.time_colour(local_now, self.HEADING_COLOUR),
        ))

        self.ax.text(
            x=self.coord[0], y=self.coord[1] + 6, s=f'{here_azimuth:.1f}°',
            transform=self.geodetic, zorder=12,
            color=self.time_colour(local_now, self.ANNOTATE_COLOUR),
        )
        self.ax.text(
            x=0.92, y=0.82,
            s=f'Qibla:\n'
            f'{distance*1e-3:,.0f} km @ {here_azimuth:.1f}°',
            transform=self.ax.transAxes, zorder=12,
            color=self.ANNOTATE_COLOUR[1],  # it's always "night" in space
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
    kaaba_hemi.plot_salah_isocurves(night=night)
    home_hemi.plot_salah_isocurves(night=night)
    home_hemi.plot_heading(utcnow=utcnow)

    return fig


def main() -> None:
    utcnow = datetime.now().astimezone(timezone.utc)
    home_coord = load_home()

    plt.style.use('dark_background')
    plot_spherical(home_coord=home_coord, utcnow=utcnow)
    plt.show()


if __name__ == '__main__':
    main()
