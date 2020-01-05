"""Statistical functions."""
from warnings import warn

import iris

import numpy as np

from ..coord_utils import UM_LATLON, UM_TIME, ensure_bounds
from ..exceptions import AeolusWarning
from ..grid import area_weights_cube
from ..subset import extract_last_year


__all__ = (
    "spatial",
    "spatial_quartiles",
    "meridional_mean",
    "minmaxdiff",
    "region_mean_diff",
    "zonal_mean",
    "last_year_mean",
)


def spatial(cube, aggr, coords=UM_LATLON):
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
    >>> spatial(my_data_cube, "mean")

    """
    ensure_bounds(cube)
    flag = all(cube.coord(c).has_bounds() for c in coords)
    aggregator = getattr(iris.analysis, aggr.upper())
    if flag and isinstance(aggregator, iris.analysis.WeightedAggregator):
        kw = {"weights": area_weights_cube(cube, normalize=True).data}
    else:
        kw = {}
    return cube.collapsed(coords, aggregator, **kw)


def spatial_quartiles(cube):
    """Calculate quartiles over horizontal coordinates."""
    warn("No weights are applied!", AeolusWarning)
    q25 = cube.collapsed(UM_LATLON, iris.analysis.PERCENTILE, percent=25)
    q75 = cube.collapsed(UM_LATLON, iris.analysis.PERCENTILE, percent=75)
    return q25, q75


def minmaxdiff(cubelist, name):
    """
    Spatial maximum minus spatial minimum for a given cube.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    name: str
        Cube name.

    Returns
    -------
    iris.cube.Cube
        Difference between the extrema with collapsed spatial dimensions.
    """
    _min = spatial(cubelist.extract_strict(name), "min")
    _max = spatial(cubelist.extract_strict(name), "max")
    diff = _max - _min
    diff.rename(f"{name}_difference")
    return diff


def region_mean_diff(cubelist, name, region_a, region_b):
    """
    Difference between averages over two regions for a given cube.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    name: str
        Cube name.
    region_a: aeolus.region.Region
        First region.
    region_b: aeolus.region.Region
        Second region.

    Returns
    -------
    iris.cube.Cube
        Difference between the region averages with collapsed spatial dimensions.
    """
    mean_a = spatial(cubelist.extract_strict(name).extract(region_a.constraint), "mean")
    mean_b = spatial(cubelist.extract_strict(name).extract(region_b.constraint), "mean")
    diff = mean_a - mean_b
    diff.rename(f"{name}_mean_diff_{region_a}_{region_b}")
    return diff


def meridional_mean(cube, lat_name=UM_LATLON[0]):
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
    coslat2d = iris.util.broadcast_to_shape(
        coslat, cube.shape, cube.coord_dims(lat_name)
    )
    cube_mean = (cube * coslat2d).collapsed(lat_name, iris.analysis.SUM) / np.sum(
        coslat
    )
    return cube_mean


def zonal_mean(cube, lon_name=UM_LATLON[1]):
    """
    Calculate cube's zonal average.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with a latitude coordinate.
    lon_name: str, optional
        Name of the longitude coordinate.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.
    """
    cube_mean = cube.collapsed(lon_name, iris.analysis.MEAN)
    return cube_mean


def last_year_mean(cube):
    """Get the time mean of over the last year."""
    last_year = extract_last_year(cube)
    return last_year.collapsed(UM_TIME, iris.analysis.MEAN)
