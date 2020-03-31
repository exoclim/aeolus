# -*- coding: utf-8 -*-
"""Functionality related to coordinates of cubes."""
from cartopy.util import add_cyclic_point

import iris
from iris.util import broadcast_to_shape, guess_coord_axis

from .exceptions import ArgumentError, NotFoundError


__all__ = (
    "UM_TIME",
    "UM_HGT",
    "UM_LEV",
    "UM_LATLON",
    "UM_Z_COORDS",
    "UM_TIME_COORDS",
    "add_cyclic_point_to_cube",
    "coord_to_cube",
    "ensure_bounds",
    "get_cube_datetimes",
    "get_dim_coord",
    "nearest_coord_value",
    "not_equal_coord_axes",
    "regrid_3d",
    "replace_z_coord",
)

UM_TIME = "time"
UM_HGT = "level_height"
UM_LATLON = ["latitude", "longitude"]
UM_SIGMA = "sigma"
UM_LEV = "model_level_number"
UM_Z_COORDS = [UM_SIGMA, UM_LEV]
UM_TIME_COORDS = ["forecast_reference_time", "forecast_period", UM_TIME]


def get_cube_datetimes(cube):
    """Get a list of `iris.cube.Cube`'s time points as `datetime.datetime`s."""
    return cube.coord("time").units.num2date(cube.coord("time").points)


def nearest_coord_value(cube, coord_name, val):
    """
    Get the nearest value of a coordinate.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with the coordinate
    coord_name: str or iris.coords.Coord
        Coordinate where to look the nearest point up
    val: int or float
        The value to find

    Returns
    -------
    int or float
        element of the coordinate array closest to the given `val`
    """
    coord = cube.coord(coord_name)
    i = coord.nearest_neighbour_index(val)
    return coord.points[i]


def coord_to_cube(cube, coord):
    """
    Convert coordinate points to a cube of the same dimension as the given cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube containing the coordinate to be broadcast.
    coord: str or iris.coords.Coord
        Coordinate to be broadcast

    Returns
    -------
    iris.cube.Cube
        Cube of broadcast coordinate
    """
    if isinstance(coord, str):
        _coord = cube.coord(coord)
    else:
        _coord = coord
    dim_map = cube.coord_dims(_coord.name())
    _data = _coord.points
    if len(dim_map) > 0:
        _data = broadcast_to_shape(_data, cube.shape, dim_map)
        dc = [(c.copy(), cube.coord_dims(c)) for c in cube.dim_coords]
        ac = [(c.copy(), cube.coord_dims(c)) for c in cube.aux_coords]
        new_cube = iris.cube.Cube(
            data=_data,
            units=_coord.units,
            long_name=_coord.name(),
            dim_coords_and_dims=dc,
            aux_coords_and_dims=ac,
        )
    else:
        new_cube = iris.cube.Cube(data=_data, standard_name=_coord.name(), units=_coord.units)
    return new_cube


def ensure_bounds(cube, coords=UM_LATLON):
    """Auto-generate bounds for cube coordinates."""
    for coord_name in coords:
        c = cube.coord(coord_name)
        if not c.has_bounds():
            if len(c.points) > 1:
                c.guess_bounds()


def not_equal_coord_axes(cube1, cube2):
    """Given 2 cubes, return axes of unequal dimensional coordinates."""
    coord_comp = iris.analysis.coord_comparison(cube1, cube2)
    neq_dim_coords = set(coord_comp["not_equal"]).intersection(set(coord_comp["dimensioned"]))
    dims = []
    for coord_pair in neq_dim_coords:
        for coord in coord_pair:
            dims.append(iris.util.guess_coord_axis(coord))
    return set(filter(None, dims))


def regrid_3d(cube, target, vert_coord=None):
    """
    Regrid a cube in the horizontal and in the vertical on to coordinates of the target cube.

    Adapted from https://github.com/LSaffin/iris-extensions

    Parameters
    ----------
    cube: iris.cube.Cube
        The cube to be regridded.
    target: iris.cube.Cube
        The cube to regrid to.
    vert_coord: str or iris.coords.Coord, optional
        The coordinate for the vertical interpolation.
        If not given, the target's z-axis `iris.coord.DimCoord` is used.

    Returns
    -------
        iris.cube.Cube
    """
    neq_axes = not_equal_coord_axes(cube, target)
    if neq_axes.intersection(["X", "Y"]):
        cube = cube.regrid(target, iris.analysis.Linear())

    # Interpolate in the vertical if needed
    if "Z" in neq_axes:
        if vert_coord is None:
            z = get_dim_coord(target, "z")
        else:
            z = target.coord(vert_coord)
        cube = cube.interpolate([(z.name(), z.points)], iris.analysis.Linear())
        ensure_bounds(cube, coords=[z])

    # Match coordinate information
    # XXX is this needed?
    # newcube = target.copy(data=cube.data)
    # newcube.rename(cube.name())
    # newcube.units = cube.units

    # Put back correct time information
    # for coord in newcube.aux_coords:
    #     if iris.util.guess_coord_axis(coord) == "T":
    #         newcube.remove_coord(coord)
    # for coord in cube.aux_coords:
    #     if iris.util.guess_coord_axis(coord) == "T":
    #         newcube.add_aux_coord(coord)
    return cube


def get_dim_coord(cube, axis):
    """
    Return a coordinate from a cube based on the axis it represents.

    Uses :py:func:`iris.util.guess_coord_axis` to heuristically match a dimensional coordinate
    with the requested axis.

    Adapted from https://github.com/LSaffin/iris-extensions

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with the desired coordinate.
    axis: str
        The co-ordinate axis to take from the cube. Must be one of X, Y, Z, T.

    Returns
    -------
    iris.coords.DimCoord
        The dimensional coordinate matching the requested axis on the given cube.

    Raises
    ------
    ArgumentError: If axis is not one of {X, Y, Z, T}.
    NotFoundError: If the cube does not contain a coord with the requested axis.
    """
    _allowed = ["X", "Y", "Z", "T"]
    axis = axis.upper()
    # If the axis supplied is not correct raise an error
    if axis not in _allowed:
        raise ArgumentError(f"Axis must be one of {_allowed}, {axis} is given.")

    # Loop over dimensional coords in the cube
    for coord in cube.dim_coords:
        # Return the coordinate if it matches the axis
        if axis == guess_coord_axis(coord):
            return coord

    # If no coordinate matches raise an error
    raise NotFoundError(f"Cube has no coordinate for axis {axis}")


def replace_z_coord(cube, promote_coord=UM_HGT, remove_coord=UM_Z_COORDS):
    """
    Replace dimensional vertical coordinate.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.
    promote_coord: str or iris.coords.Coord
        Coordinate to become a dimensional z-coordinate.
    remove_coord: list-like
        List of coordinates to remove.
        By default, model levels and sigma coordinates are removed.

    Returns
    -------
    iris.cube.Cube
        Copy of the input cube with a new vertical coordinate.
    """
    new_cube = cube.copy()
    new_cube.coord(promote_coord).bounds = None
    iris.util.promote_aux_coord_to_dim_coord(new_cube, promote_coord)
    ensure_bounds(new_cube, coords=[promote_coord])
    for coord in remove_coord:
        try:
            new_cube.remove_coord(coord)
        except iris.exceptions.CoordinateNotFoundError:
            pass
    return new_cube


def add_cyclic_point_to_cube(cube, coord=UM_LATLON[1]):
    """
    Add a cyclic point to a cube and a corresponding coordinate.

    A wrapper for `cartopy.util.add_cyclic_point()`, generalising it for iris cubes.

    Parameters
    ----------
    cube: iris.cube.Cube
        An n-dimensional cube of data to add a cyclic point to.
    coord: iris.coords.Coord or str
        A 1-dimensional coordinate which specifies the coordinate values for
        the dimension the cyclic point is to be added to. The coordinate
        values must be regularly spaced. Defaults to "longitude".

    Returns
    -------
    cyclic_cube
        The cube with a cyclic point added.
    """
    the_coord = cube.coord(coord)
    dim = cube.coord_dims(the_coord)

    cy_data, cy_coord_pnts = add_cyclic_point(cube.data, coord=the_coord.points, axis=dim[0])

    dim_coords_and_dims = [
        (coord, cube.coord_dims(coord)) for coord in cube.dim_coords if coord != the_coord
    ]
    dim_coords_and_dims.append((the_coord.copy(cy_coord_pnts), dim))

    aux_coords_and_dims = [
        (coord, cube.coord_dims(coord))
        for coord in cube.aux_coords
        if cube.coord_dims(coord) != dim
    ]
    other_kwargs = {
        key: getattr(cube, key, None)
        for key in [
            "attributes",
            "standard_name",
            "long_name",
            "var_name",
            "units",
            "cell_methods",
            "aux_factories",
            "cell_measures_and_dims",
        ]
    }
    cyclic_cube = iris.cube.Cube(
        cy_data,
        dim_coords_and_dims=dim_coords_and_dims,
        aux_coords_and_dims=aux_coords_and_dims,
        **other_kwargs,
    )
    return cyclic_cube
