import json
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

    GEODETIC_COLOUR = 'lightgreen'
    FEATURE_COLOUR = 'black'
    ANNOTATE_COLOUR = 'white'
    HEADING_COLOUR = 'red'

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

    def plot(
        self,
        night: Nightshade,
        utcnow: datetime,
    ) -> None:
        self.ax.set_title(f'{self.name} hemisphere')
        self.ax.stock_img()
        self.ax.add_feature(night, zorder=5)
        self.ax.gridlines()

        self.ax.plot(
            [self.coord[0], self.endpoint[0]],
            [self.coord[1], self.endpoint[1]],
            transform=self.geodetic, c=self.GEODETIC_COLOUR, zorder=10,
        )

        self.ax.scatter(
            [self.coord[0]], [self.coord[1]],
            transform=self.geodetic, c=self.FEATURE_COLOUR, marker='+', zorder=11,
        )

        local_now = utcnow.astimezone(self.timezone)
        self.ax.text(
            x=self.coord[0], y=self.coord[1],
            s=local_now.strftime('%Y-%m-%d %H:%M:%S %z'),
            transform=self.geodetic, rotation=270,
            color=self.ANNOTATE_COLOUR, ha='left', va='top', zorder=12,
        )

        x, y = night._geoms[0].boundary.coords.xy
        x = np.array(x)
        y = np.array(y)
        is_sunset = x < 0
        is_sunrise = x > 0

        xr, yr, zr = self.geodetic.transform_points(
            src_crs=night.crs,
            x=x, y=y,
        ).T
        self.ax.plot(
            xr[is_sunrise], yr[is_sunrise], transform=self.geodetic, c='orange', zorder=10,
        )
        self.ax.plot(
            xr[is_sunset], yr[is_sunset], transform=self.geodetic, c='purple', zorder=10,
        )

    def plot_heading(self) -> None:
        geodesic = Geodesic()  # WGS-84, not spherical - so we don't pass in class params
        (distance, home_azimuth, kabaa_azimuth), = geodesic.inverse(
            points=self.coord, endpoints=self.endpoint,
        )

        north = 0
        self.ax.plot(
            [self.coord[0], self.coord[0]],
            [self.coord[1], self.coord[1] + 20],
            transform=self.geodetic, c=self.HEADING_COLOUR, zorder=10,
        )
        self.ax.add_patch(Arc(
            xy=self.coord, width=10, height=10,
            transform=self.geodetic, color=self.HEADING_COLOUR,
            theta1=90 - home_azimuth, theta2=90 - north,
            zorder=10,
        ))

        self.ax.text(
            x=self.coord[0], y=self.coord[1] + 6, transform=self.geodetic,
            s=f'{home_azimuth:.1f}°', color=self.ANNOTATE_COLOUR, zorder=12,
        )
        self.ax.text(
            x=0.92, y=0.82, transform=self.ax.transAxes, zorder=12,
            s=f'WGS-84 elliptical geodesic\n'
            f'{distance*1e-3:,.0f} km @ {home_azimuth:.1f}°',
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
    kaaba_hemi.plot(night=night, utcnow=utcnow)
    home_hemi.plot(night=night, utcnow=utcnow)
    home_hemi.plot_heading()

    return fig


def main() -> None:
    utcnow = datetime.now().astimezone(timezone.utc)
    home_coord = load_home()

    plt.style.use('dark_background')
    plot_spherical(home_coord=home_coord, utcnow=utcnow)
    plt.show()


if __name__ == '__main__':
    main()
