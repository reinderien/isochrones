import json
import os
import shutil
import typing
from pathlib import Path
from urllib.parse import urljoin

from requests import Session

if typing.TYPE_CHECKING:
    from .types import Coord


def load_home() -> 'Coord':
    """
    Load the prayer coordinate
    """
    with open('.home.json') as f:
        coords = json.load(f)
    return coords['lon'], coords['lat']


def load_blue_marble() -> None:
    """
    Download two resolutions of the July 2004 Blue Marble image.
    File details are described in https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73751/readme.pdf
    File listing is in https://visibleearth.nasa.gov/images/73751/july-blue-marble-next-generation-w-topography-and-bathymetry/73753l
    Projection is Plate Carrée, WGS84 datum.
    Save to .resources, and make cartopy look there.
    """

    resources = Path.cwd() / '.resources'
    if not resources.is_dir():
        raise FileNotFoundError()
    os.environ['CARTOPY_USER_BACKGROUNDS'] = str(resources)

    base = 'https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73751/'

    with Session() as session:
        for filename in (
            'world.topo.bathy.200407x294x196.jpg',
            'world.topo.bathy.200407.3x5400x2700.jpg',
        ):
            file_path = resources / filename
            if not file_path.is_file():
                url = urljoin(base, filename)
                print(f'Downloading {url}...')
                with session.get(
                    url=url, stream=True,
                    headers={'Accept': 'image/jpeg'},
                ) as resp:
                    resp.raise_for_status()
                    with file_path.open(mode='wb') as f:
                        shutil.copyfileobj(fsrc=resp.raw, fdst=f)
