"""Generic calculus functions."""
import iris
from iris.analysis.calculus import differentiate


import numpy as np

from ..model import um


__all__ = (
    "d_dx",
    "d_dy",
    "d_dz",
    "deriv",
    "integrate",
)


def d_dx(cube, model=um):
    """Derivative w.r.t. x-coordinate."""
    return deriv(cube, model.x)


def d_dy(cube, model=um):
    """Derivative w.r.t. y-coordinate."""
    return deriv(cube, model.y)


def d_dz(cube, model=um):
    """Derivative w.r.t. z-coordinate."""
    return deriv(cube, model.z)


def deriv(cube, coord):
    """
    Derivative w.r.t. the given coordinate.

    Uses `iris.analysis.calculus.differentiate` and then
    interpolates the result to the grid points of the original cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube containing the given coordinate.
    coord: str or iris.coords.Coord
        Coordinate for differentiation.

    Returns
    -------
    iris.cube.Cube
        d(cube)/d(coord).

    See also
    --------
    aeolus.calc.calculus.d_dx, aeolus.calc.calculus.d_dy, aeolus.calc.calculus.d_dz
    """
    pnts = cube.coord(coord).points
    diff = differentiate(cube, coord)
    res = diff.interpolate([(coord, pnts)], iris.analysis.Linear())
    return res


def integrate(cube, coord):
    """
    Integrate the cube along a 1D coordinate using the trapezoidal rule.

    Note: `coord` must be one of the dimensional coordinates of the cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube containing the given coordinate.
    coord: str or iris.coords.Coord
        Coordinate for integration.

    Returns
    -------
    iris.cube.Cube
        Integrated cube.
    """
    # TODO: allow non-dim coordinates
    c = cube.coord(coord)
    others = [dc.name() for dc in cube.dim_coords if cube.coord_dims(dc) != cube.coord_dims(c)]
    dim = cube.coord_dims(c)[0]
    data = np.trapz(cube.data, c.points, axis=dim)
    res = next(cube.slices(others)).copy(data=data)
    res.units = cube.units * c.units
    res.remove_coord(c)
    res.rename(f"integral_of_{cube.name()}_wrt_{c.name()}")
    # ensure_bounds(cube, [c])
    # delta = iris.coords.AuxCoord(c.bounds[:, 1] - c.bounds[:, 0], units=c.units)
    # res = iris.analysis.maths.multiply(cube, delta, dim=dim).collapsed(c, iris.analysis.SUM)
    return res
