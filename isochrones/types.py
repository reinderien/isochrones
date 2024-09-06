import typing
import numpy as np

Metre = float
Radian = float
Degree = float
GeoDeg = Degree  # geographic
EclDeg = Degree  # ecliptic
Second = float

# Coordinates, unitless
Coord = tuple[float, float]

# Coordinates, degrees, unspecified whether ecliptic or geographic
CoordDeg = tuple[Degree, Degree]

CoordEclDeg = tuple[EclDeg, EclDeg]
CoordGeoDeg = tuple[GeoDeg, GeoDeg]

_FloatArray = np.ndarray[
    typing.Any, np.dtype[np.float64],
]
DegArray = _FloatArray
RadArray = _FloatArray
