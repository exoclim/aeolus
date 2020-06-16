"""Subset cubes using iris constraints."""
from datetime import timedelta

import iris

from .coord import get_cube_datetimes
from .model import um


__all__ = (
    "CM_INST_CONSTR",
    "CM_MEAN_CONSTR",
    "DimConstr",
    "extract_last_month",
    "extract_last_year",
    "l_range_constr",
)


def _select_inst(cube):
    """Return True if cube has no cell methods."""
    return len(cube.cell_methods) == 0


def _select_mean(cube):
    """Return True if cube has a "mean" cell method."""
    try:
        return cube.cell_methods[0].method == "mean"
    except IndexError:
        return False


def _dim_constr(*coords, strict=True):
    """Make an `iris.Constraint` from given dimensional coordinates."""
    coord_set = set(coords)

    def __cube_func(cube):  # noqa
        cube_dimcoords = {dc.name() for dc in cube.dim_coords}
        if strict:
            # Cube has exactly the same dim coords
            return cube_dimcoords == coord_set
        else:
            # Cube has these dim coords and possibly other
            return coord_set.issubset(cube_dimcoords)

    return iris.Constraint(cube_func=__cube_func)


def l_range_constr(h_min, h_max, units="km", coord=um.z):
    """Make a constraint on length range."""
    if units == "km":
        factor = 1e-3
    else:
        factor = 1
    return iris.Constraint(**{coord: lambda x: h_min <= (x.point * factor) <= h_max})


def extract_last_month(cube, model=um):
    """Extract time slices within the last months of a cube."""
    dt = get_cube_datetimes(cube)[-1]
    return cube.extract(
        iris.Constraint(
            **{model.t: lambda x: (x.point.year == dt.year) and (x.point.month == dt.month)}
        )
    )


def extract_last_year(cube, model=um):
    """Extract time slices within the last year of a cube."""
    dt = get_cube_datetimes(cube)[-1]
    yr_before = dt - timedelta(days=365)
    return cube.extract(iris.Constraint(**{model.t: lambda t: t.point > yr_before}))


CM_INST_CONSTR = iris.Constraint(cube_func=_select_inst)
CM_MEAN_CONSTR = iris.Constraint(cube_func=_select_mean)


class DimConstr:
    """Container for some useful dimensional constraints."""

    def __init__(self, model=um):
        """
        Initialise DimConstr.

        Parameters
        ----------
        model: aeolus.model.Model, optional
            Model class with relevant coordinate names.
        """
        self.tmyx = _dim_constr(model.t, model.l, model.y, model.x)
        self.tzyx = _dim_constr(model.t, model.z, model.y, model.x)
        self.tyx = _dim_constr(model.t, model.y, model.x)
        self.myx = _dim_constr(model.l, model.y, model.x)
        self.zyx = _dim_constr(model.z, model.y, model.x)
        self.yx = _dim_constr(model.y, model.x)
        self.yx_r = _dim_constr(model.y, model.x, strict=False)
