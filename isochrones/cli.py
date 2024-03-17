import argparse
import logging
from datetime import datetime, timezone


def utc_parse(s: str) -> datetime:
    return datetime.fromisoformat(s).astimezone(timezone.utc)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog='python -m isochrones',
        description='Islamic prayer time visualisation',
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-s', '--static', action='store_true',
        help='Static plot (animate otherwise)',
    )
    group.add_argument(
        '-f', '--time-factor', type=float, default=1,
        help='Time factor for animation; default 1. May be negative to '
             'travel back in time. Use 3600 for one hour per second.',
    )

    # Don't set the default here; take 'now' closer to time of plot
    parser.add_argument(
        '-t', '--time', type=utc_parse,
        help='ISO 8601-format timestamp of static plot or '
             'beginning of animation sequence; default now',
    )
    parser.add_argument(
        '-x', '--longitude', type=float,
        help='Home longitude in degrees (overrides .home.json)',
    )
    parser.add_argument(
        '-y', '--latitude', type=float,
        help='Home latitude in degrees (overrides .home.json)',
    )

    return parser


def main() -> None:
    args = make_parser().parse_args()

    logger = logging.getLogger(__name__)

    # This can take a bit of time, so do it locally rather than at the top
    logger.info('Loading libraries...')
    from matplotlib import pyplot as plt
    from .files import load_home, load_blue_marble
    from .plot import animate_spherical, plot_spherical

    logger.info('Loading files...')
    load_blue_marble()
    if args.longitude is None:
        home_coord = load_home()
    else:
        home_coord = args.longitude, args.latitude

    logger.info('Setting up artists...')
    if args.static:
        fig = plot_spherical(home_coord=home_coord, utcnow=args.time)
    else:
        fig, anim = animate_spherical(
            home_coord=home_coord, start_utc=args.time, time_factor=args.time_factor,
        )

    logger.info('Plotting...')
    plt.show()
