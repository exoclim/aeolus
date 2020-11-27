"""Subset cubes using iris constraints."""
from datetime import timedelta
from itertools import combinations

import iris

from .coord import get_cube_datetimes
from .model import um


__all__ = (
    "CM_INST_CONSTR",
    "CM_MEAN_CONSTR",
    "DimConstr",
    "extract_last_month",
    "extract_last_n_days",
    "l_range_constr",
    "unique_cubes",
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


def extract_last_n_days(cube, days=365, model=um):
    """Extract time slices within the last `n` days of its time dimension."""
    dt = get_cube_datetimes(cube)[-1]
    ndays_before = dt - timedelta(days=days)
    cube_sub = cube.extract(iris.Constraint(**{model.t: lambda t: t.point > ndays_before}))
    return cube_sub


CM_INST_CONSTR = iris.Constraint(cube_func=_select_inst)
CM_MEAN_CONSTR = iris.Constraint(cube_func=_select_mean)


class _ModeDimConstr:
    """Class for storing constraints in `DimConstr` grouped by their mode."""

    def __init__(self, **kwargs):
        """Assign all keyword arguments to attributes of this class."""
        self.__dict__.update(**kwargs)


class DimConstr:
    """
    Container for strict or relaxed dimensional constraints.

    Examples
    --------
    Extract cubes that have y and x dimensional coordinates (among others):
    >>> dc = DimConstr()
    >>> cubelist.extract(dc.relax.yx)

    Extract cubes that only have model levels, y and x dimensions:
    >>> dc = DimConstr()
    >>> cubelist.extract(dc.strict.myx)
    """

    def __init__(self, model=um):
        """
        Initialise DimConstr.

        Parameters
        ----------
        model: aeolus.model.Model, optional
            Model class with relevant coordinate names.
        """
        abbr_aliases = {"t": "t", "z": "z", "m": "lev", "y": "y", "x": "x"}
        for mode in ["strict", "relax"]:
            attrs = {}
            for key in ("z", "m"):
                for n in range(1, 5):
                    for seq in combinations(["t", key, "y", "x"], n):
                        model_seq = [abbr_aliases[letter] for letter in seq]
                        attrs["".join(seq)] = _dim_constr(
                            *[getattr(model, i) for i in model_seq], strict=(mode == "strict")
                        )
            setattr(self, mode, _ModeDimConstr(**attrs))


def unique_cubes(cubelist):
    """Remove duplicate cubes from `iris.cube.CubeList`."""
    if len(cubelist) > 1:
        out = iris.cube.CubeList([cubelist[0]])
        for src_cube in cubelist[1:]:
            exists = False
            for dest_cube in out:
                if src_cube == dest_cube:
                    exists = True
            if not exists:
                out.append(src_cube)
        return out
    else:
        return cubelist
