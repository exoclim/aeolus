"""Operations on geographical grid."""
from warnings import warn

import iris
from iris.analysis.cartography import wrap_lons
from iris.util import is_regular

import numpy as np

from .coord_utils import UM_LATLON, ensure_bounds
from .exceptions import AeolusWarning, LoadError


__all__ = (
    "roll_cube_0_360",
    "roll_cube_pm180",
    "area_weights_cube",
    "add_binned_lon_lat",
    "coarsen_cube",
)


def _is_longitude_global(lon_points):
    """Return True if array of longitudes covers the whole sphere."""
    dx = np.diff(lon_points)[0]  # assume regular grid
    case_0_360 = ((lon_points[0] - dx) <= 0) and ((lon_points[-1] + dx) >= 360)
    case_pm180 = ((lon_points[0] - dx) <= -180) and ((lon_points[-1] + dx) >= 180)
    return case_0_360 or case_pm180


def roll_cube_pm180(cube_in, coord_name=UM_LATLON[1], inplace=False):
    """
    Take a cube spanning 0...360 degrees in longitude and roll it to -180...180 degrees.

    Works with global model output, and in some cases for regional.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    coord_name: str, optional
        Name of the longitude coordinate.
    inplace: bool, optional
        Do this in-place or copy the cube.

    Returns
    -------
    iris.cube.Cube

    See also
    --------
    aeolus.grid.roll_cube_0_360
    """
    if inplace:
        cube = cube_in
    else:
        cube = cube_in.copy()
    xcoord = cube.coord(coord_name)
    if (xcoord.points >= 0.0).all():
        assert is_regular(xcoord), "Operation is only valid for a regularly spaced coordinate."
        if _is_longitude_global(xcoord.points):
            # Shift data symmetrically only when dealing with global cubes
            cube.data = np.roll(cube.data, len(xcoord.points) // 2, axis=-1)

        if xcoord.has_bounds():
            bounds = np.sort(wrap_lons(xcoord.bounds, -180, 360), axis=0)  # + subtract
        else:
            bounds = None
        cube.replace_coord(
            xcoord.copy(points=np.sort(wrap_lons(xcoord.points, -180, 360)), bounds=bounds)
        )
    else:
        # Nothing to do, the cube is already centered on 0 longitude
        # unless there is something wrong with longitude
        msg = f"Incorrect {coord_name} values: from {xcoord.points.min()} to {xcoord.points.max()}"
        assert ((xcoord.points >= -180.0) & (xcoord.points <= 180.0)).all(), msg
    if not inplace:
        return cube


def roll_cube_0_360(cube_in, inplace=False):
    """
    Take a cube spanning -180...180 degrees in longitude and roll it to 0...360 degrees.

    Works with global model output, and in some cases for regional.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    coord_name: str, optional
        Name of the longitude coordinate.
    inplace: bool, optional
        Do this in-place or copy the cube.

    Returns
    -------
    iris.cube.Cube

    See also
    --------
    aeolus.grid.roll_cube_pm180
    """
    if inplace:
        cube = cube_in
    else:
        cube = cube_in.copy()
    lon = cube.coord("longitude")
    if (lon.points < 0.0).any():
        add = 180
        cube.data = np.roll(cube.data, len(lon.points) // 2, axis=-1)
        if lon.has_bounds():
            bounds = lon.bounds + add
        else:
            bounds = None
        cube.replace_coord(lon.copy(points=lon.points + add, bounds=bounds))
    if not inplace:
        return cube


def area_weights_cube(cube, r_planet=None, normalize=False):
    """
    Create a cube of area weights for an arbitrary planet.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates
    r_planet: float, optional
        Radius of the planet.
    normalize: bool, optional
        Normalize areas.

    Returns
    -------
    iris.cube.Cube
        Cube of area weights with the same metadata as the input cube
    """
    cube = cube.copy()
    ensure_bounds(cube)
    aw = iris.analysis.cartography.area_weights(cube, normalize=normalize)
    if normalize:
        aw = cube.copy(data=aw)
        aw.rename("normalized_grid_cell_area")
        aw.units = "1"
    else:
        cs = cube.coord_system("CoordSystem")
        if r_planet is not None:
            r = r_planet
        elif cs is not None:
            r = cs.semi_major_axis
        else:
            try:
                r = cube.attributes["planet_conf"].radius
            except (KeyError, LoadError):
                warn("Using default Earth radius", AeolusWarning)
                r = None

        if r is not None:
            aw *= (r / iris.fileformats.pp.EARTH_RADIUS) ** 2
        aw = cube.copy(data=aw)
        aw.rename("grid_cell_area")
        aw.units = "m**2"
    return aw


def _cell_bounds(points, bound_position=0.5):
    """
    Calculate coordinate cell boundaries.

    Taken from SciTools iris package.

    Parameters
    ----------
    points: numpy.array
        One-dimensional array of uniformy spaced values of shape (M,)
    bound_position: bool, optional
        The desired position of the bounds relative to the position
        of the points.

    Returns
    -------
    bounds: numpy.array
        Array of shape (M+1,)

    Examples
    --------
    >>> a = np.arange(-1, 2.5, 0.5)
    >>> a
    array([-1. , -0.5,  0. ,  0.5,  1. ,  1.5,  2. ])
    >>> cell_bounds(a)
    array([-1.25, -0.75, -0.25,  0.25,  0.75,  1.25,  1.75,  2.25])

    See Also
    --------
    aeolus.grid._cell_centres
    """
    assert points.ndim == 1, "Only 1D points are allowed"
    diffs = np.diff(points)
    if not np.allclose(diffs, diffs[0]):
        warn("_cell_bounds() is supposed to work only for uniformly spaced points", AeolusWarning)
    delta = diffs[0] * bound_position
    bounds = np.concatenate([[points[0] - delta], points + delta])
    return bounds


def _cell_centres(bounds, bound_position=0.5):
    """
    Calculate coordinate cell centres.

    Taken from SciTools iris package.

    Parameters
    ----------
    bounds: numpy.array
        One-dimensional array of cell boundaries of shape (M,)
    bound_position: bool, optional
        The desired position of the bounds relative to the position
        of the points.

    Returns
    -------
    centres: numpy.array
        Array of shape (M-1,)

    Examples
    --------
    >>> a = np.arange(-1, 3., 1.)
    >>> a
    array([-1,  0,  1,  2])
    >>> cell_centres(a)
    array([-0.5,  0.5,  1.5])

    See Also
    --------
    aeolus.grid._cell_bounds
    """
    assert bounds.ndim == 1, "Only 1D points are allowed"
    deltas = np.diff(bounds) * bound_position
    centres = bounds[:-1] + deltas
    return centres


def add_binned_lon_lat(cube, lon_bins, lat_bins, coord_names=UM_LATLON, inplace=False):
    """
    Add binned longitude and latitude as auxiliary coordinates to a cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    lon_bins: array-like
        Longitude bins.
    lat_bins: array-like
        Latitude bins.
    coord_names: list, optional
        List of latitude and longitude labels.
    inplace: bool, optional
        Do this in-place or copy the cube.

    Returns
    -------
    iris.cube.Cube
    """
    if inplace:
        cube_out = cube
    else:
        cube_out = cube.copy()
    for name, target_points in zip(coord_names, (lat_bins, lon_bins)):
        binned_points = np.digitize(cube_out.coord(name).points, target_points)
        binned_points = np.clip(binned_points, 0, len(target_points) - 1)
        new_coord = iris.coords.AuxCoord(binned_points, long_name=f"{name}_binned")
        cube_out.add_aux_coord(new_coord, cube_out.coord_dims(name))
    return cube_out


def coarsen_cube(cube, lon_bins, lat_bins, coord_names=UM_LATLON, inplace=False):
    """
    Block-average cube in longitude and latitude.

    Note: no weighting is applied!

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    lon_bins: array-like
        Longitude bins.
    lat_bins: array-like
        Latitude bins.
    coord_names: list, optional
        List of latitude and longitude labels.
    inplace: bool, optional
        Do this in-place or copy the cube.

    Returns
    -------
    iris.cube.Cube
    """
    if inplace:
        cube_out = cube
    else:
        cube_out = cube.copy()
    add_binned_lon_lat(cube_out, lon_bins, lat_bins, coord_names=coord_names, inplace=True)

    # To avoid oversampling on the edges, extract subset within the boundaries of target coords
    for coord, target_points in zip(coord_names, (lat_bins, lon_bins)):
        cube_out = cube_out.extract(
            iris.Constraint(**{coord: lambda p: target_points.min() <= p <= target_points.max()})
        )

    for coord in coord_names:
        cube_out = cube_out.aggregated_by([f"{coord}_binned"], iris.analysis.MEAN)

    for coord, target_points in zip(coord_names, (lat_bins, lon_bins)):
        dim = cube_out.coord_dims(coord)
        units = cube_out.coord(coord).units
        cube_out.remove_coord(coord)
        aux = cube_out.coord(f"{coord}_binned")
        new_points = target_points[aux.points]
        #         if len(aux.points) < len(target_points):
        #             new_points =
        new_coord = iris.coords.DimCoord.from_coord(aux.copy(points=new_points, bounds=None))
        cube_out.remove_coord(f"{coord}_binned")
        new_coord.rename(coord)
        new_coord.units = units
        cube_out.add_dim_coord(new_coord, dim)

    return cube_out
