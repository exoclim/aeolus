# -*- coding: utf-8 -*-
"""Subset cubes using iris constraints."""
from datetime import timedelta
from itertools import combinations

from iris import Constraint
from iris.cube import CubeList

from .coord import get_cube_datetimes
from .model import um


__all__ = (
    "CM_INST_CONSTR",
    "CM_MEAN_CONSTR",
    "DimConstr",
    "extract_after_n_days",
    "extract_between_days",
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

    return Constraint(cube_func=__cube_func)


def l_range_constr(h_min, h_max, units="km", coord=um.z):
    """Make a constraint on length range."""
    if units == "km":
        factor = 1e-3
    else:
        factor = 1
    return Constraint(**{coord: lambda x: h_min <= (x.point * factor) <= h_max})


def extract_last_month(cube, model=um):
    """Extract time slices within the last months of a cube."""
    dt = get_cube_datetimes(cube)[-1]
    return cube.extract(
        Constraint(**{model.t: lambda x: (x.point.year == dt.year) and (x.point.month == dt.month)})
    )


def extract_last_n_days(cube, days=365, model=um):
    """Extract time slices within the last `n` days of its time dimension."""
    dt = get_cube_datetimes(cube)[-1]
    ndays_before = dt - timedelta(days=days)
    cube_sub = cube.extract(Constraint(**{model.t: lambda t: t.point > ndays_before}))
    return cube_sub


def extract_between_days(cube, day_start, day_end, model=um):
    """Extract a cube subset between `day_start` and `day_end` of its time coordinate."""
    dt_s = get_cube_datetimes(cube)[0]
    lbound = dt_s + timedelta(days=day_start)
    ubound = dt_s + timedelta(days=day_end)
    cube_sub = cube.extract(Constraint(**{model.t: lambda t: lbound <= t.point <= ubound}))
    return cube_sub


def extract_after_n_days(cube, days=365, model=um):
    """Extract time slices after the given number of `days`."""
    dt = get_cube_datetimes(cube)[0]
    ndays_after = dt + timedelta(days=days)
    cube_sub = cube.extract(Constraint(**{model.t: lambda t: t.point > ndays_after}))
    return cube_sub


CM_INST_CONSTR = Constraint(cube_func=_select_inst)
CM_MEAN_CONSTR = Constraint(cube_func=_select_mean)


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
        abbr_aliases = {"t": "t", "z": "z", "m": "lev", "p": "p", "y": "y", "x": "x"}
        for mode in ["strict", "relax"]:
            attrs = {}
            for key in ("z", "m", "p"):
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
        out = CubeList([cubelist[0]])
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


class CellMethodConstraint(Constraint):
    """Provides a simple Cube-CellMethod based :class:`Constraint`."""

    _names = ("method", "coord_names", "intervals", "comments")

    def __init__(self, method, coord_names=None, intervals=None, comments=None):
        """
        Initialise a CellMethod-based constraint.

        Provides a Cube cell method based :class:`Constraint`, which matches
        against attributes of each of the cell methods in the cube, if there
        are any. See :class:`iris.coords.CellMethod` for more info.

        The constraint will succeed if *any* of the provided attributes match.

        Example usage::

            iris.CellMethodConstraint(method='mean')

            iris.CellMethodConstraint(method='mean',
                coord_names=("mesh_coordinates",))

        .. note:: Cell method names are case sensitive.

        """
        self.method = method
        self.coord_names = coord_names
        self.intervals = intervals
        self.comments = comments
        # self.cell_method = iris.coords.CellMethod(
        #     method, coords=coord_names, intervals=intervals, comments=comments
        # )
        super().__init__(cube_func=self._cube_func)

    def __eq__(self, other):
        eq = type(other) is CellMethodConstraint and all(
            getattr(self, attname) == getattr(other, attname) for attname in self._names
        )
        return eq

    def __hash__(self):
        # Must re-define if you overload __eq__ : Use object identity.
        return id(self)

    def _cube_func(self, cube):
        match = False
        try:
            cell_methods = cube.cell_methods
            for cell_method in cell_methods:
                for name in self._names:
                    expected = getattr(self, name)
                    actual = getattr(cell_method, name)
                    # Currently relaxed / permissive condition
                    if actual == expected:
                        match = True
                        break
        except AttributeError:
            match = False
        return match

    def __repr__(self):
        names = []
        for name in self._names:
            value = getattr(self, name)
            if value != "none":
                names.append("{}={!r}".format(name, value))
        return "{}({})".format(self.__class__.__name__, ", ".join(names))
