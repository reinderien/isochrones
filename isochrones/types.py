import typing
import numpy as np

Metre = float
Radian = float
Degree = float
GeoDeg = Degree  # geographic
EclDeg = Degree  # ecliptic
Second = float

Coord = tuple[float, float]
CoordDeg = tuple[Degree, Degree]
CoordEclipticDeg = tuple[EclDeg, EclDeg]
CoordGeoDeg = tuple[GeoDeg, GeoDeg]

DegArray = np.ndarray[
    typing.Any, np.dtype[Degree],
]
RadArray = np.ndarray[
    typing.Any, np.dtype[Radian],
]
