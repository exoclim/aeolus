"""Operations on geographical grid."""
import iris

import numpy as np


def roll_cube_e2w(cube_in, inplace=False):
    """
    Take a cube which goes longitude 0-360 back to -180-180.

    Works with global model output, and in some cases for regional.
    """
    if inplace:
        cube = cube_in
    else:
        cube = cube_in.copy()
    lon = cube.coord("longitude")
    if (lon.points >= 0.0).all():
        if (lon.points <= 360.0).all():
            subtract = -180.0
            cube.data = np.roll(cube.data, len(lon.points) // 2, axis=-1)
        else:
            # because for regional runs, UM output
            # can have longitudes > 360
            # In this case no data roll needed
            subtract = -360.0

        if lon.has_bounds():
            bounds = lon.bounds + subtract
        else:
            bounds = None
        cube.replace_coord(lon.copy(points=lon.points + subtract, bounds=bounds))
    else:
        # Nothing to do, the cube is already centered on 0 longitude
        # unless there is something wrong with longitude
        msg = "Incorrect longitude values:" f"from {lon.points.min()} to {lon.points.max()}"
        assert ((lon.points >= -180.0) & (lon.points <= 180.0)).all(), msg
    if not inplace:
        return cube


def area_weights_cube(cube, r_planet=None):
    """
    Create a cube of area weights for an arbitrary planet.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates
    r_planet: float
        Radius of planet

    Returns
    -------
    iris.cube.Cube
        Cube of area weights with the same metadata as the input cube
    """
    cube = cube.copy()
    for dim_coord in cube.dim_coords:
        if not dim_coord.has_bounds():
            dim_coord.guess_bounds()
    aw = iris.analysis.cartography.area_weights(cube)
    if r_planet is not None:
        aw *= (r_planet / iris.fileformats.pp.EARTH_RADIUS) ** 2  # FIXME
    aw = cube.copy(data=aw)
    aw.units = "m**2"
    return aw
