"""Operations on geographical grid."""
from warnings import warn

import iris
from iris.analysis.cartography import wrap_lons
from iris.util import is_regular

import numpy as np

from .exceptions import AeolusWarning


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
