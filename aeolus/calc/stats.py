"""Statistical functions."""
from warnings import warn

import iris
from iris.util import broadcast_to_shape

import numpy as np

from .calculus import integrate
from ..coord import UM_HGT, UM_LATLON, UM_TIME, area_weights_cube, ensure_bounds
from ..exceptions import AeolusWarning
from ..subset import extract_last_year


__all__ = (
    "spatial",
    "spatial_quartiles",
    "meridional_mean",
    "minmaxdiff",
    "region_mean_diff",
    "zonal_mean",
    "last_year_mean",
    "vertical_cumsum",
    "vertical_mean",
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
    coslat2d = iris.util.broadcast_to_shape(coslat, cube.shape, cube.coord_dims(lat_name))
    cube_mean = (cube * coslat2d).collapsed(lat_name, iris.analysis.SUM) / np.sum(coslat)
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


def vertical_mean(cube, coord=UM_HGT, weight_by=None):
    """
    Vertical mean of a cube with optional weighting.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube to average.
    coord: str or iris.coords.Coord
        Coordinate of the dimension over which the averaging is done.
    weight_by: str or iris.coords.Coord or iris.cube.Cube, optional
        Coordinate of the given cube or another cube used for weighting.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.
    """
    if weight_by is None:
        vmean = cube.collapsed(coord, iris.analysis.MEAN)
    else:
        if isinstance(weight_by, (str, iris.coords.Coord)):
            vmean = cube.collapsed(
                UM_HGT, iris.analysis.MEAN, weights=cube.coord(weight_by).points.squeeze()
            )
        elif isinstance(weight_by, iris.cube.Cube):
            a_copy = cube.copy()
            b_copy = weight_by.copy()
            a_copy.coord(coord).bounds = None
            b_copy.coord(coord).bounds = None
            prod = b_copy * a_copy
            vmean = integrate(prod, coord) / integrate(weight_by, coord)
        else:
            raise ValueError(f"unrecognised type of weight_by: {type(weight_by)}")
    vmean.rename(f"vertical_mean_of_{cube.name()}")
    return vmean


def vertical_cumsum(cube, coord=UM_HGT):
    """
    Vertical cumulative sum of a cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.
    coord: str or iris.coords.Coord
        Coordinate of the dimension over which the cumulative sum is calculated.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.
    """
    c = cube.coord(coord).copy()
    dim = cube.coord_dims(c)
    c.guess_bounds()
    z_weights = broadcast_to_shape(c.bounds[:, 1] - c.bounds[:, 0], cube.shape, dim)
    data = np.nancumsum(cube.data * z_weights, axis=dim[0])
    res = cube.copy(data=data)
    res.rename(f"vertical_cumulative_sum_of_{cube.name()}")
    res.units = cube.units * c.units
    return res
