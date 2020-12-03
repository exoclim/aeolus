"""Generic calculus functions."""
from cf_units import Unit

import iris
from iris.analysis.calculus import _coord_cos, _curl_differentiate, _curl_regrid, differentiate

import numpy as np

from .meta import update_metadata
from ..const import get_planet_radius
from ..exceptions import NotYetImplementedError
from ..model import um


__all__ = (
    "d_dx",
    "d_dy",
    "d_dz",
    "deriv",
    "div_h",
    "integrate",
)


def d_dx(cube, model=um):
    """Calculate a derivative w.r.t. x-coordinate."""
    return deriv(cube, model.x)


def d_dy(cube, model=um):
    """Calculate a derivative w.r.t. y-coordinate."""
    return deriv(cube, model.y)


def d_dz(cube, model=um):
    """Calculate a derivative w.r.t. z-coordinate."""
    return deriv(cube, model.z)


def deriv(cube, coord):
    """
    Calculate a derivative w.r.t. the given coordinate.

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


@update_metadata(name="horizontal_divergence")
def div_h(i_cube, j_cube, r_planet=None, model=um):
    r"""
    Calculate horizontal divergence.

    Note: currently works only in spherical coordinates.

    Parameters
    ----------
    i_cube:
        i-component (e.g. u-wind)
    j_cube:
        j-component (e.g. v-wind)
    r_planet: float, optional
        Radius of the planet (m). If not given, an attempt is made
        to get it from the cube metadata.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of horizontal divergence.

    Notes
    -----
    Divergence in spherical coordinates is defined as

    .. math::

        \nabla\cdot \vec A = \frac{1}{r cos \theta} (
        \frac{\partial \vec A_\lambda}{\partial \lambda}
        + \frac{\partial}{\partial \theta}
        (\vec A_\theta cos \theta))

    where \lambda is longitude, \theta is latitude.
    """
    x_coord = i_cube.coord(model.x)
    y_coord = i_cube.coord(model.y)

    y_dim = i_cube.coord_dims(y_coord)[0]

    horiz_cs = i_cube.coord_system("CoordSystem")

    # Check for spherical coords
    spherical_coords = isinstance(
        horiz_cs, (iris.coord_systems.GeogCS, iris.coord_systems.RotatedGeogCS)
    )
    if spherical_coords:
        # Get the radius of the planet
        if r_planet is None:
            r = get_planet_radius(i_cube)
        else:
            r = r_planet
        r_unit = Unit("m")

        lon_coord = x_coord.copy()
        lat_coord = y_coord.copy()
        lon_coord.convert_units("radians")
        lat_coord.convert_units("radians")
        lat_cos_coord = _coord_cos(lat_coord)

        # j-component: \frac{\partial}{\partial \theta} (\vec A_\theta cos \theta))
        temp = iris.analysis.maths.multiply(j_cube, lat_cos_coord, y_dim)
        djcos_dtheta = _curl_differentiate(temp, lat_coord)
        prototype_diff = djcos_dtheta

        # i-component: \frac{\partial \vec A_\lambda}{\partial \lambda}
        d_i_cube_dlambda = _curl_differentiate(i_cube, lon_coord)
        d_i_cube_dlambda = _curl_regrid(d_i_cube_dlambda, prototype_diff)
        new_lat_coord = d_i_cube_dlambda.coord(model.y)
        new_lat_cos_coord = _coord_cos(new_lat_coord)
        lat_dim = d_i_cube_dlambda.coord_dims(new_lat_coord)[0]

        # Sum and divide
        div = iris.analysis.maths.divide(
            iris.analysis.maths.add(d_i_cube_dlambda, djcos_dtheta),
            r * new_lat_cos_coord,
            dim=lat_dim,
        )
        div.units /= r_unit
        div = div.regrid(i_cube, iris.analysis.Linear())
    else:
        raise NotYetImplementedError("Non-spherical coordinates are not implemented yet.")
    return div


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
