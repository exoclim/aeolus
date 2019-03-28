"""Operations on geographical grid."""
import iris
from iris.analysis.cartography import wrap_lons

import numpy as np


def _is_longitude_global(lon_points):
    """Return True if array of longitudes covers the whole sphere."""
    dx = np.diff(lon_points)[0]  # assume regular grid
    case_0_360 = ((lon_points[0] - dx) <= 0) and ((lon_points[-1] + dx) >= 360)
    case_pm180 = ((lon_points[0] - dx) <= -180) and ((lon_points[-1] + dx) >= 180)
    return case_0_360 or case_pm180


def roll_cube_e2w(cube_in, coord_name="longitude", inplace=False):
    """
    Take a cube which goes longitude 0-360 back to -180-180.

    Works with global model output, and in some cases for regional.
    """
    if inplace:
        cube = cube_in
    else:
        cube = cube_in.copy()
    xcoord = cube.coord(coord_name)
    if (xcoord.points >= 0.0).all():
        if _is_longitude_global(xcoord.points):
            # Shift data symmetrically only when dealing with global cubes
            cube.data = np.roll(cube.data, len(xcoord.points) // 2, axis=-1)

        if xcoord.has_bounds():
            bounds = wrap_lons(xcoord.bounds, -180, 360)  # + subtract
        else:
            bounds = None
        cube.replace_coord(xcoord.copy(points=wrap_lons(xcoord.points, -180, 360), bounds=bounds))
    else:
        # Nothing to do, the cube is already centered on 0 longitude
        # unless there is something wrong with longitude
        msg = f"Incorrect {coord_name} values: from {xcoord.points.min()} to {xcoord.points.max()}"
        assert ((xcoord.points >= -180.0) & (xcoord.points <= 180.0)).all(), msg
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
    aw.rename("grid_cell_area")
    aw.units = "m**2"
    return aw
