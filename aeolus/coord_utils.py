# -*- coding: utf-8 -*-
"""Functionality related to coordinates of cubes."""
import iris
from iris.util import broadcast_to_shape


def get_cube_datetimes(cube):
    """Get a list of `iris.cube.Cube`'s time points as `datetime.datetime`s."""
    return cube.coord("time").units.num2date(cube.coord("time").points)


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
    new_cube: iris.cube.Cube
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
            standard_name=_coord.name(),
            dim_coords_and_dims=dc,
            aux_coords_and_dims=ac,
        )
    else:
        new_cube = iris.cube.Cube(data=_data, standard_name=_coord.name(), units=_coord.units)
    return new_cube
