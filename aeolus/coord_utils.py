# -*- coding: utf-8 -*-
"""Functionality related to coordinates of cubes."""
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
    "get_cube_datetimes",
    "nearest_coord_value",
    "coord_to_cube",
    "ensure_bounds",
    "not_equal_coord_axes",
    "regrid_3d",
    "get_dim_coord",
)

UM_TIME = "time"
UM_LEV = "model_level_number"
UM_HGT = "level_height"
UM_LATLON = ["latitude", "longitude"]
UM_Z_COORDS = ["sigma", UM_LEV]
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


# def z_interp_cube(cube, z=None, z_coord_name=UM_HGT, replace_z_coord=True):
#     """
#     Interpolate cube to z points in vertical and use height as dim coord.
#
#     Parameters
#     ----------
#     cube: iris.cube.Cube
#         Input cube
#     z: numpy.array, optional
#         Array of z-points as target for interpolation.
#     z_coord_name: str, optional
#         Vertical coordinate for interpolation.
#     replace_z_coord: bool, optional
#         Replace model levels with level height.
#
#     Returns
#     -------
#     iris.cube.Cube
#         Interpolated cube with updated coordinate data
#     """
#     cube_out = cube.copy()
#     # Remove model level numbers and other redundant coords
#     if replace_z_coord:
#         for coord in UM_Z_COORDS:
#             try:
#                 cube_out.remove_coord(coord)
#             except iris.exceptions.CoordinateNotFoundError:
#                 pass
#         # and use level height as dim coord
#         try:
#             iris.util.promote_aux_coord_to_dim_coord(cube_out, z_coord_name)
#         except ValueError:
#             cube_out.coord(z_coord_name).bounds = None
#             cube_out.coord(z_coord_name).guess_bounds()
#             iris.util.promote_aux_coord_to_dim_coord(cube_out, z_coord_name)
#
#     if z is not None:
#         height_target = [(z_coord_name, z)]
#         cube_out = cube_out.interpolate(height_target, iris.analysis.Linear())
#
#     # Code below is for rescaling z-coordinate
#     # lh = cube_out.coord(z_coord_name)
#     # dim = cube_out.coord_dims(z_coord_name)[0]
#     # cube_out.remove_coord(z_coord_name)
#     # lh = lh.copy(lh.points * zscale, bounds=None)
#     # cube_out.add_dim_coord(lh, dim)
#     return cube_out
