"""Operations on geographical grid."""
import numpy as np


def roll_cube_e2w(cube_in, inplace=False):
    """
    Takes a cube which goes longitude 0-360 back to -180-180.
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
