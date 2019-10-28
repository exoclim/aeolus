"""Functions to calculate commonly used variables in atmospheric science."""
import iris

import numpy as np

from .coord_utils import ensure_bounds, UM_LATLON
from .grid import area_weights_cube


def calc_spatial(cube, aggr, coords=UM_LATLON):
    """
    Calculate spatial statistic with geographic grid weights.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    aggr: str
        Statistical aggregator (see iris.analysis for available aggregators).
    coords: list, optional
        List of names of spatial coordinates.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.

    Examples
    --------
    >>> calc_spatial(my_data_cube, "mean")

    """
    ensure_bounds(cube)
    flag = all([cube.coord(c).has_bounds() for c in coords])
    aggregator = getattr(iris.analysis, aggr.upper())
    if flag and isinstance(aggregator, iris.analysis.WeightedAggregator):
        kw = dict(weights=area_weights_cube(cube, normalize=True).data)
    else:
        kw = {}
    return cube.collapsed(coords, aggregator, **kw)


def calc_meridional_mean(cube, lat_name=UM_LATLON[0]):
    """
    Calculate cube's meridional average.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with a latitude coordinate.
    lat_name: str, optional
        Name of the latitude coordinate.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.
    """
    coslat = np.cos(np.deg2rad(cube.coord(lat_name).points))
    coslat2d = iris.util.broadcast_to_shape(coslat, cube.shape, (0,))
    cube_mean = (cube * coslat2d).collapsed(lat_name, iris.analysis.SUM) / np.sum(coslat)
    return cube_mean
