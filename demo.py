import json
from typing import NamedTuple

import cartopy.crs
import cartopy.feature.nightshade
import cartopy.mpl.geoaxes
import matplotlib.patches
from matplotlib import pyplot as plt

KAABA_COORD = (39.826167, 21.4225)


def load_home() -> tuple[float, float]:
    with open('.home.json') as f:
        coords = json.load(f)
    return coords['lon'], coords['lat']


class Hemisphere(NamedTuple):
    name: str
    coord: tuple[float, float]
    crs: cartopy.crs.Orthographic
    ax: cartopy.mpl.geoaxes.GeoAxes

    GEODETIC_COLOUR = 'green'
    FEATURE_COLOUR = 'black'

    @classmethod
    def make(
        cls,
        name: str,
        index: int,
        coord: tuple[float, float],
        figure: plt.Figure,
        n_axes: int = 2,
        globe: cartopy.crs.Globe | None = None,
    ) -> 'Hemisphere':
        crs = cartopy.crs.Orthographic(
            central_longitude=coord[0],
            central_latitude=coord[1],
            globe=globe,
        )
        ax = figure.add_subplot(1, n_axes, index, projection=crs)
        return cls(name=name, coord=coord, crs=crs, ax=ax)

    def plot(
        self,
        endpoint: tuple[float, float],
        geodetic: cartopy.crs.Geodetic,
        night: cartopy.feature.nightshade.Nightshade,
    ) -> None:
        self.ax.set_title(f'{self.name} hemisphere')

        self.ax.stock_img()
        self.ax.plot(
            [self.coord[0], endpoint[0]],
            [self.coord[1], endpoint[1]],
            transform=geodetic, c=self.GEODETIC_COLOUR,
        )
        self.ax.add_feature(night)
        self.ax.gridlines()
        self.ax.scatter(
            [self.coord[0]], [self.coord[1]],
            transform=geodetic, c=self.FEATURE_COLOUR, marker='+', zorder=10,
        )


class HeadingHemisphere(Hemisphere):
    HEADING_COLOUR = 'red'

    def plot_heading(
        self,
        geodetic: cartopy.crs.Geodetic,
        geodesic_azimuth: float,
    ) -> None:
        north = 90
        self.ax.plot(
            [self.coord[0], self.coord[0]],
            [self.coord[1], self.coord[1] + 20],
            transform=geodetic, c=self.HEADING_COLOUR,
        )
        self.ax.add_patch(matplotlib.patches.Arc(
            xy=self.coord, transform=geodetic, color=self.HEADING_COLOUR,
            width=10, height=10, theta1=geodesic_azimuth, theta2=north,
        ))


def plot(
    home_coord: tuple[float, float],
    geodesic_azimuth: float,
) -> plt.Figure:
    fig: plt.Figure = plt.figure()

    home_hemi = HeadingHemisphere.make(
        name='Home', index=1, coord=home_coord, figure=fig)
    kaaba_hemi = Hemisphere.make(
        name='Kaaba', index=2, coord=KAABA_COORD, figure=fig, globe=home_hemi.crs.globe)

    geodetic = cartopy.crs.Geodetic(globe=home_hemi.crs.globe)
    night = cartopy.feature.nightshade.Nightshade()

    kaaba_hemi.plot(endpoint=home_hemi.coord, geodetic=geodetic, night=night)
    home_hemi.plot(endpoint=kaaba_hemi.coord, geodetic=geodetic, night=night)
    home_hemi.plot_heading(geodetic=geodetic, geodesic_azimuth=geodesic_azimuth)

    return fig


def main() -> None:
    home_coord = load_home()

    plot(
        home_coord=home_coord,
        geodesic_azimuth=10,  # wrong
    )

    plt.show()


if __name__ == '__main__':
    main()
