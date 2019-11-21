"""Subset cubes using iris constraints."""
from datetime import timedelta

import iris

from .coord_utils import UM_HGT, UM_LATLON, UM_LEV, UM_TIME, get_cube_datetimes


__all__ = (
    "CM_MEAN_CONSTR",
    "DIM_CONSTR_TMYX",
    "DIM_CONSTR_TZYX",
    "DIM_CONSTR_TYX",
    "DIM_CONSTR_MYX",
    "DIM_CONSTR_ZYX",
    "DIM_CONSTR_YX",
    "DIM_CONSTR_YX_R",
    "extract_last_month",
    "extract_last_year",
    "l_range_constr",
)


def _select_mean(cube):
    """Return True if cube has a "mean" cell method."""
    try:
        return cube.cell_methods[0].method == "mean"
    except IndexError:
        return False


def _dim_constr(*coords, strict=True):
    """Make an `iris.Constraint` from given dimensional coordinates."""
    # make a cube function with local variables
    def __cube_func(cube):  # noqa
        cube_dimcoords = {dc.name() for dc in cube.dim_coords}
        if strict:
            # Cube has exactly the same dim coords
            return cube_dimcoords == set(coords)
        else:
            # Cube has these dim coords and possibly other
            return cube_dimcoords.issubset(coords)

    return iris.Constraint(cube_func=__cube_func)


def l_range_constr(h_min, h_max, units="km", coord=UM_HGT):
    """Make a constraint on length range."""
    if units == "km":
        factor = 1e-3
    else:
        factor = 1
    return iris.Constraint(**{coord: lambda x: h_min <= (x.point * factor) <= h_max})


def extract_last_month(cube):
    """Extract time slices within the last months of a cube."""
    dt = get_cube_datetimes(cube)[-1]
    return cube.extract(
        iris.Constraint(
            **{UM_TIME: lambda x: (x.point.year == dt.year) and (x.point.month == dt.month)}
        )
    )


def extract_last_year(cube):
    """Extract time slices within the last year of a cube."""
    dt = get_cube_datetimes(cube)[-1]
    yr_before = dt - timedelta(days=365)
    return cube.extract(iris.Constraint(**{UM_TIME: lambda t: t.point > yr_before}))


CM_MEAN_CONSTR = iris.Constraint(cube_func=_select_mean)

DIM_CONSTR_TMYX = _dim_constr(UM_TIME, UM_LEV, *UM_LATLON)
DIM_CONSTR_TZYX = _dim_constr(UM_TIME, UM_HGT, *UM_LATLON)
DIM_CONSTR_TYX = _dim_constr(UM_TIME, *UM_LATLON)
DIM_CONSTR_MYX = _dim_constr(UM_LEV, *UM_LATLON)
DIM_CONSTR_ZYX = _dim_constr(UM_HGT, *UM_LATLON)
DIM_CONSTR_YX = _dim_constr(*UM_LATLON)
DIM_CONSTR_YX_R = _dim_constr(*UM_LATLON, strict=False)
