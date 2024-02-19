import json
from datetime import datetime, tzinfo
from typing import NamedTuple

import pytz
from cartopy.crs import Geodetic, Orthographic
from cartopy.feature.nightshade import Nightshade
from cartopy.geodesic import Geodesic
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib import pyplot as plt
from matplotlib.patches import Arc

KAABA_COORD = (39.826167, 21.4225)
KAABA_TIMEZONE = pytz.timezone('Asia/Riyadh')


def load_home() -> tuple[float, float]:
    with open('.home.json') as f:
        coords = json.load(f)
    return coords['lon'], coords['lat']


class Hemisphere(NamedTuple):
    name: str
    coord: tuple[float, float]
    crs: Orthographic
    geodetic: Geodetic
    ax: GeoAxes
    timezone: tzinfo

    GEODETIC_COLOUR = 'lightgreen'
    FEATURE_COLOUR = 'black'
    ANNOTATE_COLOUR = 'white'

    @classmethod
    def make(
        cls,
        name: str,
        index: int,
        coord: tuple[float, float],
        figure: plt.Figure,
        timezone: tzinfo,
        n_axes: int = 2,
        geodetic: Geodetic | None = None,
    ) -> 'Hemisphere':
        crs = Orthographic(
            central_longitude=coord[0],
            central_latitude=coord[1], globe=geodetic and geodetic.globe,
        )
        ax = figure.add_subplot(1, n_axes, index, projection=crs)
        return cls(
            name=name, coord=coord, crs=crs, ax=ax, timezone=timezone,
            geodetic=geodetic or Geodetic(globe=crs.globe))

    def plot(
        self,
        endpoint: tuple[float, float],
        night: Nightshade,
        home_now: datetime,
    ) -> None:
        self.ax.set_title(f'{self.name} hemisphere')
        self.ax.stock_img()
        self.ax.add_feature(night, zorder=5)
        self.ax.gridlines()

        self.ax.plot(
            [self.coord[0], endpoint[0]],
            [self.coord[1], endpoint[1]],
            transform=self.geodetic, c=self.GEODETIC_COLOUR, zorder=10,
        )

        self.ax.scatter(
            [self.coord[0]], [self.coord[1]],
            transform=self.geodetic, c=self.FEATURE_COLOUR, marker='+', zorder=11,
        )

        local_now = home_now.astimezone(self.timezone)
        self.ax.text(
            x=self.coord[0], y=self.coord[1],
            s=local_now.strftime('%Y-%m-%d %H:%M:%S %z'),
            transform=self.geodetic, rotation=270,
            color=self.ANNOTATE_COLOUR, ha='left', va='top', zorder=12,
        )


class HeadingHemisphere(Hemisphere):
    HEADING_COLOUR = 'red'

    def plot_heading(
        self,
        geodesic_azimuth: float,
        geodesic_distance: float,
    ) -> None:
        north = 0
        self.ax.plot(
            [self.coord[0], self.coord[0]],
            [self.coord[1], self.coord[1] + 20],
            transform=self.geodetic, c=self.HEADING_COLOUR, zorder=10,
        )
        self.ax.add_patch(Arc(
            xy=self.coord, width=10, height=10,
            transform=self.geodetic, color=self.HEADING_COLOUR,
            theta1=90 - geodesic_azimuth, theta2=90 - north,
            zorder=10,
        ))

        self.ax.text(
            x=self.coord[0], y=self.coord[1] + 6, transform=self.geodetic,
            s=f'{geodesic_azimuth:.1f}°', color=self.ANNOTATE_COLOUR, zorder=12,
        )
        self.ax.text(
            x=0.92, y=0.82, transform=self.ax.transAxes, zorder=12,
            s=f'WGS-84 elliptical geodesic\n'
            f'{geodesic_distance*1e-3:,.0f} km @ {geodesic_azimuth:.1f}°',
        )


def plot_spherical(
    home_coord: tuple[float, float],
    geodesic_azimuth: float,
    geodesic_distance: float,
    home_now: datetime,
) -> plt.Figure:
    fig: plt.Figure = plt.figure()

    home_hemi = HeadingHemisphere.make(
        name='Home', index=1, coord=home_coord, figure=fig, timezone=home_now.tzinfo)
    kaaba_hemi = Hemisphere.make(
        name='Kaaba', index=2, coord=KAABA_COORD, figure=fig, timezone=KAABA_TIMEZONE,
        geodetic=home_hemi.geodetic,
    )

    night = Nightshade()

    kaaba_hemi.plot(endpoint=home_hemi.coord, night=night, home_now=home_now)
    home_hemi.plot(endpoint=kaaba_hemi.coord, night=night, home_now=home_now)
    home_hemi.plot_heading(
        geodesic_azimuth=geodesic_azimuth,
        geodesic_distance=geodesic_distance,
    )

    return fig


def main() -> None:
    home_now = datetime.now().astimezone()
    home_coord = load_home()

    geodesic = Geodesic()  # WGS-84
    (distance, home_azimuth, kabaa_azimuth), = geodesic.inverse(
        points=home_coord, endpoints=KAABA_COORD,
    )

    plot_spherical(
        geodesic_azimuth=home_azimuth,
        geodesic_distance=distance, home_coord=home_coord,
        home_now=home_now,
    )

    plt.show()


if __name__ == '__main__':
    main()
