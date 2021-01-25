"""Statistical functions."""
from warnings import warn

import iris
from iris.util import broadcast_to_shape

import numpy as np

from .calculus import integrate
from ..coord import area_weights_cube, ensure_bounds
from ..exceptions import AeolusWarning
from ..model import um
from ..subset import extract_last_n_days


__all__ = (
    "cumsum",
    "last_n_day_mean",
    "meridional_mean",
    "minmaxdiff",
    "normalize_cube",
    "region_mean_diff",
    "spatial",
    "spatial_mean",
    "spatial_quartiles",
    "time_mean",
    "vertical_mean",
    "zonal_mean",
)


def cumsum(cube, axis, axis_weights=False, model=um):
    """
    Cumulative sum of a cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.
    axis: str
        Coordinate axis of operation (t|z|y|x).
    axis_weights: bool, optional
        If True, multiply data by the coordinate spacings.
    model: aeolus.model.Model, optional
        Model class with a relevant coordinate name.

    Returns
    -------
    iris.cube.Cube
        Cube of cumulative sums with the same dimensions as the input cube.
    """
    try:
        c = cube.coord(getattr(model, axis)).copy()
    except (AttributeError, iris.exceptions.CoordinateNotFoundError):
        c = cube.coord(axis=axis).copy()
    dim = cube.coord_dims(c)
    if axis_weights:
        if not c.has_bounds():
            c.guess_bounds()
        weights = broadcast_to_shape(c.bounds[:, 1] - c.bounds[:, 0], cube.shape, dim)
        data = cube.data * weights
        units = cube.units * c.units
    else:
        data = cube.data
        units = cube.units
    data = np.nancumsum(data, axis=dim[0])
    res = cube.copy(data=data)
    res.rename(f"cumulative_sum_of_{cube.name()}_along_{axis}")
    res.units = units
    return res


def last_n_day_mean(cube, days=365, model=um):
    """Average the cube over the last `n` days of its time dimension."""
    cube_sub = time_mean(extract_last_n_days(cube, days=days, model=model), model=model)
    return cube_sub


def meridional_mean(cube, model=um):
    """
    Calculate cube's meridional average.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with a latitude coordinate.
    model: aeolus.model.Model, optional
        Model class with a relevant coordinate name.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.
    """
    lat_name = model.y
    coslat = np.cos(np.deg2rad(cube.coord(lat_name).points))
    coslat2d = iris.util.broadcast_to_shape(coslat, cube.shape, cube.coord_dims(lat_name))
    cube_mean = (cube * coslat2d).collapsed(lat_name, iris.analysis.SUM) / np.sum(coslat)
    return cube_mean


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
    _min = spatial(cubelist.extract_cube(name), "min")
    _max = spatial(cubelist.extract_cube(name), "max")
    diff = _max - _min
    diff.rename(f"{name}_difference")
    return diff


def normalize_cube(cube):
    """
    Normalize cube data, i.e. make the values range from 0 to 1.

    .. math::
        z_i = (x_i - min(x)) / (max(x) - min(x))

    Parameters
    ----------
    cube: iris.cube.Cube
        The input cube.
    """
    cube_min = cube.collapsed(cube.dim_coords, iris.analysis.MIN)
    cube_max = cube.collapsed(cube.dim_coords, iris.analysis.MAX)
    norm = (cube - cube_min) / (cube_max - cube_min)
    norm.rename(f"normalized_{cube.name()}")
    return norm


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
    mean_a = spatial_mean(cubelist.extract_cube(name).extract(region_a.constraint))
    mean_b = spatial_mean(cubelist.extract_cube(name).extract(region_b.constraint))
    diff = mean_a - mean_b
    diff.rename(f"{name}_mean_diff_{region_a}_{region_b}")
    return diff


def spatial(cube, aggr, model=um):
    """
    Calculate spatial statistic with geographic grid weights.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    aggr: str
        Statistical aggregator (see iris.analysis for available aggregators).
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.

    Examples
    --------
    >>> spatial(my_data_cube, "mean")

    """
    ensure_bounds(cube, coords=("x", "y"), model=model)
    coords = (model.y, model.x)
    flag = all(cube.coord(c).has_bounds() for c in coords)
    aggregator = getattr(iris.analysis, aggr.upper())
    if flag and isinstance(aggregator, iris.analysis.WeightedAggregator):
        kw = {"weights": area_weights_cube(cube, normalize=True).data}
    else:
        kw = {}
    return cube.collapsed(coords, aggregator, **kw)


def spatial_mean(cube, model=um):
    """Shortcut for spatial(cube, "mean")."""
    return spatial(cube, "mean", model=model)


def spatial_quartiles(cube, model=um):
    """Calculate quartiles over horizontal coordinates."""
    warn("No weights are applied!", AeolusWarning)
    q25 = cube.collapsed((model.y, model.x), iris.analysis.PERCENTILE, percent=25)
    q75 = cube.collapsed((model.y, model.x), iris.analysis.PERCENTILE, percent=75)
    return q25, q75


def time_mean(cube, model=um):
    """Time average of a cube."""
    return cube.collapsed(model.t, iris.analysis.MEAN)


def vertical_mean(cube, weight_by=None, model=um):
    """
    Vertical mean of a cube with optional weighting.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube to average.
    weight_by: str or iris.coords.Coord or iris.cube.Cube, optional
        Coordinate of the given cube or another cube used for weighting.
    model: aeolus.model.Model, optional
        Model class with a relevant coordinate name.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.
    """
    coord = model.z
    if weight_by is None:
        vmean = cube.collapsed(coord, iris.analysis.MEAN)
    else:
        if isinstance(weight_by, (str, iris.coords.Coord)):
            weights = broadcast_to_shape(
                cube.coord(weight_by).points.squeeze(), cube.shape, cube.coord_dims(weight_by)
            )
            vmean = cube.collapsed(coord, iris.analysis.MEAN, weights=weights)
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


def zonal_mean(cube, model=um):
    """
    Calculate cube's zonal average.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with a latitude coordinate.
    model: aeolus.model.Model, optional
        Model class with a relevant coordinate name.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.
    """
    lon_name = model.x
    cube_mean = cube.collapsed(lon_name, iris.analysis.MEAN)
    return cube_mean
