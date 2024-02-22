import json


def load_home() -> tuple[float, float]:
    """
    Load the prayer coordinate
    """
    with open('.home.json') as f:
        coords = json.load(f)
    return coords['lon'], coords['lat']
