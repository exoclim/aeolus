"""Operations in tidally-locked coordinates."""
import iris
from iris.analysis.cartography import _meshgrid, rotate_pole, rotate_winds
from iris.coord_systems import RotatedGeogCS

import numpy as np

from ..coord import get_xy_coords
from ..exceptions import BadCoordinateError
from ..model import um

__all__ = (
    "regrid_to_rotated_pole_coordinates",
    "regrid_to_tidally_locked_coordinates",
    "rotate_winds_to_tidally_locked_coordinates",
)


def regrid_to_rotated_pole_coordinates(cube, pole_lon, pole_lat, model=um):
    """
    Regrid a cube to rotated-pole coordinates.

    Note: non-dimensional coordinates are lost.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.
    pole_lon: float
        New North Pole longitude.
    pole_lat: float
        New North Pole latitude.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate and variable names.
    """
    _cs = cube.coord_system()
    xcoord, ycoord = get_xy_coords(cube, model=model)
    lon2d, lat2d = _meshgrid(xcoord.points, ycoord.points)
    r_lons, r_lats = rotate_pole(lon2d, lat2d, pole_lon, pole_lat)

    r_lon_coord = iris.coords.AuxCoord(
        r_lons.flatten(),
        units=xcoord.units,
        standard_name=xcoord.standard_name,
        var_name=xcoord.var_name,
        coord_system=_cs,
    )
    r_lat_coord = iris.coords.AuxCoord(
        r_lats.flatten(),
        units=ycoord.units,
        standard_name=ycoord.standard_name,
        var_name=ycoord.var_name,
        coord_system=_cs,
    )
    non_xy_coords = [(i, cube.coord_dims(i)) for i in cube.dim_coords if i not in [xcoord, ycoord]]
    # assume latitude and longitude are the rightmost dimensions
    lat_dim = len(non_xy_coords)
    if lat_dim != cube.coord_dims(ycoord)[0]:
        raise BadCoordinateError("Ensure latitude and longitude are the rightmost dimensions.")
    cube_flat = iris.cube.Cube(
        cube.core_data().reshape((*cube.shape[:lat_dim], np.product(cube.shape[lat_dim:]))),
        aux_coords_and_dims=[
            *non_xy_coords,
            (r_lat_coord, lat_dim),
            (r_lon_coord, lat_dim),
        ],
        units=cube.units,
    )
    for coord, _ in non_xy_coords:
        iris.util.promote_aux_coord_to_dim_coord(cube_flat, coord)

    out = cube_flat.regrid(cube, iris.analysis.UnstructuredNearest())
    return out


def regrid_to_tidally_locked_coordinates(cube, pole_lon=0, pole_lat=0, model=um):
    """
    Regrid a cube to tidally locked coordinates.

    The substellar and antistellar point become the north and south pole, respectively.
    By default, the substellar point is assumed to be at (0, 0).

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.
    pole_lon: float, optional
        New North Pole longitude.
    pole_lat: float, optional
        New North Pole latitude.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate and variable names.
    """
    return regrid_to_rotated_pole_coordinates(
        cube, pole_lon=pole_lon, pole_lat=pole_lat, model=model
    )


def rotate_winds_to_tidally_locked_coordinates(u, v, pole_lon=0, pole_lat=0):
    """
    Rotate the horizontal wind components to tidally locked coordinates.

    The substellar and antistellar point become the north and south pole, respectively.
    By default, the substellar point is assumed to be at (0, 0).

    Parameters
    ----------
    u: iris.cube.Cube
        An instance of :class:`iris.cube.Cube` that contains the x-component of the vector.
    v: iris.cube.Cube
        An instance of :class:`iris.cube.Cube` that contains the y-component of the vector.
    pole_lon: float, optional
        New North Pole longitude.
    pole_lat: float, optional
        New North Pole latitude.
    """
    tl_cs = RotatedGeogCS(pole_lon, pole_lat, ellipsoid=u.coord_system())
    tl_u, tl_v = rotate_winds(u, v, tl_cs)
    return tl_u, tl_v
