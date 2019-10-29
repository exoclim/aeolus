# -*- coding: utf-8 -*-
"""Functionality related to coordinates of cubes."""
import iris
from iris.util import broadcast_to_shape


__all__ = (
    "get_cube_datetimes",
    "nearest_coord_value",
    "coord_to_cube",
    "z_interp_cube",
    "ensure_bounds",
)

UM_TIME_COORDS = ["forecast_reference_time", "forecast_period", "time"]
UM_Z_COORDS = ["sigma", "model_level_number"]
UM_HEIGHT = "level_height"
UM_LATLON = ["latitude", "longitude"]


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


def z_interp_cube(cube, z=None, z_coord_name=UM_HEIGHT, replace_z_coord=True):
    """
    Interpolate cube to z points in vertical and use height as dim coord.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube
    z: numpy.array, optional
        Array of z-points as target for interpolation.
    z_coord_name: str, optional
        Vertical coordinate for interpolation.
    replace_z_coord: bool, optional
        Replace model levels with level height.

    Returns
    -------
    iris.cube.Cube
        Interpolated cube with updated coordinate data
    """
    cube_out = cube.copy()
    # Remove model level numbers and other redundant coords
    if replace_z_coord:
        for coord in UM_Z_COORDS:
            try:
                cube_out.remove_coord(coord)
            except iris.exceptions.CoordinateNotFoundError:
                pass
        # and use level height as dim coord
        try:
            iris.util.promote_aux_coord_to_dim_coord(cube_out, z_coord_name)
        except ValueError:
            cube_out.coord(z_coord_name).bounds = None
            cube_out.coord(z_coord_name).guess_bounds()
            iris.util.promote_aux_coord_to_dim_coord(cube_out, z_coord_name)

    if z is not None:
        height_target = [(z_coord_name, z)]
        cube_out = cube_out.interpolate(height_target, iris.analysis.Linear())

    # lh = cube_out.coord(z_coord_name)
    # dim = cube_out.coord_dims(z_coord_name)[0]
    # cube_out.remove_coord(z_coord_name)
    # lh = lh.copy(lh.points * zscale, bounds=None)
    # cube_out.add_dim_coord(lh, dim)
    return cube_out
