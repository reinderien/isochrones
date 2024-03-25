import typing
import numpy as np

FloatArray = np.ndarray[typing.Any, np.dtype[np.float64]]

Metre = float
Degree = float
GeoDeg = Degree
EclDeg = Degree

Coord = tuple[float, float]
CoordDeg = tuple[Degree, Degree]
CoordEclipticDeg = tuple[EclDeg, EclDeg]
CoordGeoDeg = tuple[GeoDeg, GeoDeg]
