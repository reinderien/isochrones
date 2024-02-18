import json

import cartopy
from cartopy.feature.nightshade import Nightshade
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib import pyplot as plt
from matplotlib.patches import Arc

KAABA = (39.826167, 21.4225)


def load_home() -> tuple[float, float]:
    with open('.home.json') as f:
        coords = json.load(f)
    return coords['lon'], coords['lat']


def main() -> None:
    home = load_home()

    crs_home = cartopy.crs.Orthographic(
        central_longitude=home[0],
        central_latitude=home[1],
    )
    sphere_globe = crs_home.globe
    geodetic = cartopy.crs.Geodetic(globe=sphere_globe)

    crs_kabaa = cartopy.crs.Orthographic(
        central_longitude=KAABA[0],
        central_latitude=KAABA[1],
        globe=sphere_globe,
    )

    # higher-resolution geodesics - not needed for orthographic since geodesics are straight
    # crs_home.threshold *= 0.1
    # crs_kabaa.threshold *= 0.1

    geodesic_azimuth = 10  # bad

    fig: plt.Figure = plt.figure()
    axes: list[GeoAxes] = [
        fig.add_subplot(1, 2, i, projection=crs)
        for i, crs in enumerate((crs_home, crs_kabaa), start=1)
    ]
    ax_home, ax_kabaa = axes
    ax_home.set_title('Home hemisphere')
    ax_kabaa.set_title('Kabaa hemisphere')

    heading_colour = 'red'
    geodetic_colour = 'green'
    feature_colour = 'black'

    night = Nightshade()

    for ax, origin in zip(axes, (home, KAABA)):
        ax.stock_img()
        ax.plot(
            [home[0], KAABA[0]], [home[1], KAABA[1]], transform=geodetic, c=geodetic_colour,
        )
        ax.add_feature(night)
        ax.gridlines()
        ax.scatter(
            [origin[0]], [origin[1]], transform=geodetic, c=feature_colour, marker='+', zorder=10,
        )

    north = 90
    ax_home.plot(
        [home[0], home[0]],
        [home[1], home[1] + 20],
        transform=geodetic, c=heading_colour,
    )
    ax_home.add_patch(Arc(
        xy=home, transform=geodetic, color=heading_colour,
        width=10, height=10, theta1=geodesic_azimuth, theta2=north,
    ))

    plt.show()


if __name__ == '__main__':
    main()
