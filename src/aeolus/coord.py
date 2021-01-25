# -*- coding: utf-8 -*-
"""Functionality related to coordinates of cubes."""
from datetime import timedelta
from warnings import warn

from cartopy.util import add_cyclic_point

import iris
from iris.analysis.cartography import wrap_lons
from iris.coord_categorisation import _months_in_season, add_categorised_coord
from iris.exceptions import CoordinateNotFoundError as CoNotFound
from iris.experimental import stratify
from iris.util import broadcast_to_shape, guess_coord_axis, is_regular

import numpy as np

from .const import get_planet_radius
from .exceptions import AeolusWarning, ArgumentError, BadCoordinateError, NotFoundError
from .model import um


__all__ = (
    "CoordContainer",
    "add_binned_coord",
    "add_cyclic_point_to_cube",
    "add_planet_calendar",
    "area_weights_cube",
    "check_coords",
    "coarsen_cube",
    "coord_delta_to_cube",
    "coord_to_cube",
    "ensure_bounds",
    "get_cube_datetimes",
    "get_cube_rel_days",
    "get_dim_coord",
    "get_xy_coords",
    "isel",
    "interp_all_to_pres_lev",
    "interp_to_cube_time",
    "interp_to_pres_lev",
    "interp_to_single_pres_lev",
    "nearest_coord_value",
    "not_equal_coord_axes",
    "regrid_3d",
    "replace_z_coord",
    "roll_cube_0_360",
    "roll_cube_pm180",
    "vertical_cross_section_area",
    "volume_weights_cube",
)


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
    aeolus.coord._cell_centres
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
    aeolus.coord._cell_bounds
    """
    assert bounds.ndim == 1, "Only 1D points are allowed"
    deltas = np.diff(bounds) * bound_position
    centres = bounds[:-1] + deltas
    return centres


def _is_longitude_global(lon_points):
    """Return True if array of longitudes covers the whole sphere."""
    dx = np.diff(lon_points)[0]  # assume regular grid
    case_0_360 = ((lon_points[0] - dx) <= 0) and ((lon_points[-1] + dx) >= 360)
    case_pm180 = ((lon_points[0] - dx) <= -180) and ((lon_points[-1] + dx) >= 180)
    return case_0_360 or case_pm180


class CoordContainer:
    """
    Coordinate container.

    Attributes
    ----------
    {x,y,z,t}: iris.coord.Coord
        Coordinates in the respective dimensions
    """

    def __init__(self, cubes, model=um):
        """
        Instantiate an `AtmosFlow` object.

        Parameters
        ----------
        cubes: iris.cube.CubeList
            Atmospheric fields with necessary coordinates.
        model: aeolus.model.Model, optional
            Model class with relevant coordinate and variable names.
        """
        check_coords(cubes)
        # self.x = cubes[0].coord(model.x)
        self.model = model
        for axis in ["x", "y", "z", "t"]:
            try:
                setattr(self, axis, cubes[0].coord(axis=axis, dim_coords=True))
            except CoNotFound:
                pass


def add_binned_coord(cube, coord_name, bins):
    """
    Bin coordinate points and add them as an auxiliary coordinate to a cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with the given coordinate.
    coord_name: str or iris.coords.Coord
        Coordinate name.
    bins: array-like
        Coordinate bins.

    Returns
    -------
    iris.cube.Cube
    """
    cube_out = cube.copy()
    binned_points = np.digitize(cube_out.coord(coord_name).points, bins)
    binned_points = np.clip(binned_points, 0, len(bins) - 1)
    new_coord = iris.coords.AuxCoord(binned_points, long_name=f"{coord_name}_binned")
    cube_out.add_aux_coord(new_coord, cube_out.coord_dims(coord_name))
    return cube_out


def add_cyclic_point_to_cube(cube, coord=um.x):
    """
    Add a cyclic point to a cube and a corresponding coordinate.

    A wrapper for `cartopy.util.add_cyclic_point()`, generalising it for iris cubes.

    Parameters
    ----------
    cube: iris.cube.Cube
        An n-dimensional cube of data to add a cyclic point to.
    coord: iris.coords.Coord or str
        A 1-dimensional coordinate which specifies the coordinate values for
        the dimension the cyclic point is to be added to. The coordinate
        values must be regularly spaced. Defaults to the "x"-coordinate.

    Returns
    -------
    cyclic_cube
        The cube with a cyclic point added.
    """
    the_coord = cube.coord(coord)
    dim = cube.coord_dims(the_coord)

    # TODO: use core_data()?
    cy_data, cy_coord_pnts = add_cyclic_point(cube.data, coord=the_coord.points, axis=dim[0])

    dim_coords_and_dims = [
        (coord, cube.coord_dims(coord)) for coord in cube.dim_coords if coord != the_coord
    ]
    dim_coords_and_dims.append((the_coord.copy(cy_coord_pnts), dim))

    aux_coords_and_dims = [
        (coord, cube.coord_dims(coord))
        for coord in cube.aux_coords
        if cube.coord_dims(coord) != dim
    ]
    other_kwargs = {
        key: getattr(cube, key, None)
        for key in [
            "attributes",
            "standard_name",
            "long_name",
            "var_name",
            "units",
            "cell_methods",
            "aux_factories",
            "cell_measures_and_dims",
        ]
    }
    cyclic_cube = iris.cube.Cube(
        cy_data,
        dim_coords_and_dims=dim_coords_and_dims,
        aux_coords_and_dims=aux_coords_and_dims,
        **other_kwargs,
    )
    return cyclic_cube


def add_planet_calendar(
    cube,
    time_coord,
    days_in_year,
    days_in_month,
    days_in_day,
    run_start_day=0,
    seasons=("djf", "mam", "jja", "son"),
    planet="planet",
):
    """
    Add an auxiliary time axis with the non-Earth period lengths.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.
    time_coord: iris.coords.Coord or str
        Original time coordinate of the cube.
    days_in_year: int or float
        Number of Earth days in one year on the given planet.
    days_in_month: int or float
        Number of Earth days in one month on the given planet.
    days_in_day: int or float
        Number of Earth days in one day on the given planet (e.g. ~16 for Titan).
    run_start_day: int or float, optional
        Earth day of the start of the simulation.
    seasons: tuple, optional
        Sequences of letters corresponding to month names.
    planet: str, optional
        Name of the planet to be used to name the new coordinate.
    """

    def rel_day(coord, value):
        """Get the relative number of the day."""
        start = coord.units.num2date(coord.points[0])
        current = coord.units.num2date(value)
        iday = run_start_day + (current - start).days
        return iday

    def determine_season(coord, value):
        """Determine season from the month number."""
        assert coord.name() == f"{planet}_month"
        for season in seasons:
            if value + 1 in _months_in_season(season):
                return season

    new_coords = {
        "year": lambda c, v: rel_day(c, v) // days_in_year,
        "month": lambda c, v: (rel_day(c, v) % days_in_year) // days_in_month,
        "day": lambda c, v: (rel_day(c, v) % days_in_month) // days_in_day,
        "season": determine_season,
    }
    for key, op in new_coords.items():
        if key == "season":
            coord = f"{planet}_month"
        else:
            coord = time_coord
        add_categorised_coord(cube, f"{planet}_{key}", coord, op)


def area_weights_cube(cube, r_planet=None, normalize=False, model=um):
    """
    Create a cube of area weights for an arbitrary planet.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates
    r_planet: float, optional
        Radius of the planet (m). If not given, an attempt is made
        to get it from the cube metadata.
    normalize: bool, optional
        Normalize areas.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.

    Returns
    -------
    iris.cube.Cube
        Cube of area weights with the same metadata as the input cube
    """
    cube = cube.copy()
    ensure_bounds(cube, model=model)
    aw = iris.analysis.cartography.area_weights(cube, normalize=normalize)
    if normalize:
        aw = cube.copy(data=aw)
        aw.rename("normalized_grid_cell_area")
        aw.units = "1"
    else:
        if r_planet is None:
            r = get_planet_radius(cube)
        else:
            r = r_planet
        aw *= (r / iris.fileformats.pp.EARTH_RADIUS) ** 2
        aw = cube.copy(data=aw)
        aw.rename("grid_cell_area")
        aw.units = "m**2"
    return aw


def check_coords(cubes):
    """Check the cubes coordinates for consistency."""
    # get the names of all coords binned into useful comparison groups
    coord_comparison = iris.analysis._dimensional_metadata_comparison(*cubes)

    bad_coords = coord_comparison["ungroupable_and_dimensioned"]
    if bad_coords:
        raise BadCoordinateError(
            "Coordinates found in one cube that describe "
            "a data dimension which weren't in the other "
            "cube ({}), try removing this coordinate.".format(
                ", ".join(group.name() for group in bad_coords)
            )
        )

    bad_coords = coord_comparison["resamplable"]
    if bad_coords:
        raise BadCoordinateError(
            "Some coordinates are different ({}), "
            "consider resampling.".format(", ".join(group.name() for group in bad_coords))
        )


def coarsen_cube(cube, lon_bins, lat_bins, model=um):
    """
    Block-average cube in longitude and latitude.

    Note: no weighting is applied!

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    lon_bins: array-like
        Longitude bins.
    lat_bins: array-like
        Latitude bins.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.

    Returns
    -------
    iris.cube.Cube
    """
    cube_out = cube.copy()
    coord_names = [model.y, model.x]
    cube_out = add_binned_coord(cube_out, model.x, lon_bins)
    cube_out = add_binned_coord(cube_out, model.y, lat_bins)

    # To avoid oversampling on the edges, extract subset within the boundaries of target coords
    for coord, target_points in zip(coord_names, (lat_bins, lon_bins)):
        cube_out = cube_out.extract(
            iris.Constraint(**{coord: lambda p: target_points.min() <= p <= target_points.max()})
        )

    for coord in coord_names:
        cube_out = cube_out.aggregated_by([f"{coord}_binned"], iris.analysis.MEAN)

    for coord, target_points in zip(coord_names, (lat_bins, lon_bins)):
        dim = cube_out.coord_dims(coord)
        units = cube_out.coord(coord).units
        cube_out.remove_coord(coord)
        aux = cube_out.coord(f"{coord}_binned")
        new_points = target_points[aux.points]
        new_coord = iris.coords.DimCoord.from_coord(aux.copy(points=new_points, bounds=None))
        cube_out.remove_coord(f"{coord}_binned")
        new_coord.rename(coord)
        new_coord.units = units
        cube_out.add_dim_coord(new_coord, dim)

    return cube_out


def coord_delta_to_cube(cube, coord, normalize=False):
    """
    Convert 1D coordinate spacings to a cube of the same dimension as the given cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube containing the coordinate to be broadcast.
    coord: str or iris.coords.Coord
        Coordinate to be broadcast.
    normalize: bool, optional
        Normalize the data.

    Returns
    -------
    iris.cube.Cube
        Cube of broadcast coordinate deltas.
    """
    # Extract the coordinate and its dimension number
    if isinstance(coord, str):
        _coord = cube.coord(coord).copy()
    else:
        _coord = coord.copy()
    dim_map = cube.coord_dims(_coord.name())
    # Ensure the coordinate has bounds
    if not _coord.has_bounds():
        _coord.guess_bounds()

    _data = _coord.bounds[:, 1] - _coord.bounds[:, 0]
    if normalize:
        _data_max = np.nanmax(np.abs(_data))
        _data /= _data_max
        units = "1"
        prefix = "normalized_"
    else:
        units = _coord.units
        prefix = ""
    if len(dim_map) > 0:
        _data = broadcast_to_shape(_data, cube.shape, dim_map)
        dc = [(c.copy(), cube.coord_dims(c)) for c in cube.dim_coords]
        ac = [(c.copy(), cube.coord_dims(c)) for c in cube.aux_coords]
        new_cube = iris.cube.Cube(
            data=_data,
            units=units,
            dim_coords_and_dims=dc,
            aux_coords_and_dims=ac,
        )
    else:
        new_cube = iris.cube.Cube(data=_data, units=units)
    new_cube.rename(f"{prefix}delta_{_coord.name()}")
    return new_cube


def coord_to_cube(cube, coord):
    """
    Convert coordinate points to a cube of the same dimension as the given cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube containing the coordinate to be broadcast.
    coord: str or iris.coords.Coord
        Coordinate to be broadcast.

    Returns
    -------
    iris.cube.Cube
        Cube of broadcast coordinate.
    """
    if isinstance(coord, str):
        _coord = cube.coord(coord)
    else:
        _coord = coord
    dim_map = cube.coord_dims(_coord.name())
    _data = _coord.points
    if len(dim_map) > 0:
        _data = broadcast_to_shape(_data, cube.shape, dim_map)
        dc = [(c.copy(), cube.coord_dims(c)) for c in cube.dim_coords]
        ac = [(c.copy(), cube.coord_dims(c)) for c in cube.aux_coords]
        new_cube = iris.cube.Cube(
            data=_data,
            units=_coord.units,
            long_name=_coord.name(),
            dim_coords_and_dims=dc,
            aux_coords_and_dims=ac,
        )
    else:
        new_cube = iris.cube.Cube(data=_data, standard_name=_coord.name(), units=_coord.units)
    return new_cube


def ensure_bounds(cube, coords=("x", "y"), model=um):
    """Auto-generate bounds for cube coordinates."""
    for coord in coords:
        try:
            c_name = getattr(model, coord)
        except (AttributeError, TypeError):
            c_name = coord
        c = cube.coord(c_name)
        if not c.has_bounds():
            if len(c.points) > 1:
                c.guess_bounds()


def get_cube_datetimes(cube, model=um):
    """
    Convert points of a cube's time coordinate to datetimes.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube containing a time coordinate.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.

    Returns
    -------
    numpy.array
        Array of datetime-like objects.
    """
    return cube.coord(model.t).units.num2date(cube.coord(model.t).points)


def get_cube_rel_days(cube, model=um):
    """
    Convert points of a cube's time coordinate to relative number of days.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube containing a time coordinate.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.

    Returns
    -------
    numpy.array
        Array of relative days.
    """
    dts = get_cube_datetimes(cube, model=model)
    days = ((dts - dts[0]) / timedelta(days=1)).astype(np.float64)
    return days


def get_dim_coord(cube, axis):
    """
    Return a coordinate from a cube based on the axis it represents.

    Uses :py:func:`iris.util.guess_coord_axis` to heuristically match a dimensional coordinate
    with the requested axis.

    Adapted from https://github.com/LSaffin/iris-extensions

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with the desired coordinate.
    axis: str
        The co-ordinate axis to take from the cube. Must be one of X, Y, Z, T.

    Returns
    -------
    iris.coords.DimCoord
        The dimensional coordinate matching the requested axis on the given cube.

    Raises
    ------
    ArgumentError: If axis is not one of {X, Y, Z, T}.
    NotFoundError: If the cube does not contain a coord with the requested axis.
    """
    _allowed = ["X", "Y", "Z", "T"]
    axis = axis.upper()
    # If the axis supplied is not correct raise an error
    if axis not in _allowed:
        raise ArgumentError(f"Axis must be one of {_allowed}, {axis} is given.")

    # Loop over dimensional coords in the cube
    for coord in cube.dim_coords:
        # Return the coordinate if it matches the axis
        if axis == guess_coord_axis(coord):
            return coord

    # If no coordinate matches raise an error
    raise NotFoundError(f"Cube has no coordinate for axis {axis}")


def get_xy_coords(cube, model=um):
    """
    Return X and Y coordinate tuple for a given `cube` using names from a given `model` container.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate and variable names.
    """
    y_coord = cube.coord(model.y)
    x_coord = cube.coord(model.x)
    return (x_coord, y_coord)


def interp_all_to_pres_lev(cubelist, levels, interpolator=None, model=um):
    """
    Interpolate all cubes within a cubelist to the given set of pressure levels.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        List of cubes, including a cube of pressure.
    levels: array-like
        Sequence of pressure levels (same units as the units of pressure cube in `cubelist`).
    interpolator: callable or None
        The interpolator to use when computing the interpolation. See relevel() docs for more.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.CubeList
        List of cubes interpolated to pressure level(s).
    """
    pres = cubelist.extract_cube(model.pres)
    cl_out = iris.cube.CubeList()
    for cube in cubelist:
        if cube != pres:
            cube_plev = interp_to_pres_lev(
                cubelist, cube.name(), levels, interpolator=interpolator, model=model
            )
            cube_plev.coord(model.pres).attributes = {}
            cl_out.append(cube_plev)
    return cl_out


def interp_to_cube_time(cube_src, cube_tgt, model=um):
    """
    Linearly interpolate `cube_src` to `cube_tgt` along the time dimension (`model.t`).

    Forecast period is copied from `cube_tgt`.

    Parameters
    ----------
    cube_src: iris.cube.Cube
        Cube to interpolate.
    cube_tgt: iris.cube.Cube
        Cube with time dimension to interpolate to.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube with the time dimension equal to `cube_tgt`.
    """
    target = [(model.t, cube_tgt.coord(model.t).points)]
    out = cube_src.interpolate(target, iris.analysis.Linear())
    # Replace time coordinates, because interpolation removes bounds.
    for coord in ["t", "fcst_prd", "fcst_ref"]:
        try:
            out.replace_coord(cube_tgt.coord(getattr(model, coord)))
        except (AttributeError, CoNotFound):
            pass
    return out


def interp_to_pres_lev(cubelist, constraint, levels, interpolator=None, model=um):
    """
    Interpolate a cube to pressure level(s) using stratify relevel() function.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        List of cubes containing a cube extractable by `constraint` and a cube of pressure.
    constraint: str or iris.Constraint
        Variable name or constraint to extract a cube from `cubelist`.
    levels: array-like
        Sequence of pressure levels (same units as the units of pressure cube in `cubelist`).
    interpolator: callable or None
        The interpolator to use when computing the interpolation. See relevel() docs for more.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of `varname` interpolated to pressure level(s).
    """
    cube = cubelist.extract_cube(constraint)
    pres = cubelist.extract_cube(model.pres)
    cube_plev = stratify.relevel(cube, pres, levels, axis=model.z, interpolator=interpolator)
    cube_plev.coord(model.pres).attributes = {}
    return iris.util.squeeze(cube_plev)


def interp_to_single_pres_lev(
    cubelist, constraint, p_ref_frac=0.5, const=None, interpolator=None, model=um
):
    """
    Interpolate the field defined by `constraint` to a single pressure level.

    The level is found from the given fraction of the reference surface pressure.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input cubelist.
    constraint: str or iris.Constraint
        Variable name or constraint to extract a cube from `cubelist`.
    p_ref_frac: float, optional
        Fraction of reference surface pressure at which the estimate is made.
        The default value is 0.1, which for an Earth-like atmosphere means 100 hPa.
    const: aeolus.const.const.ConstContainer, optional
        If not given, constants are attempted to be retrieved from
        attributes of a cube in the cube list.
    interpolator: callable or None
        The interpolator to use when computing the interpolation. See relevel() docs for more.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube on a single pressure level.
    """
    if const is None:
        const = cubelist[0].attributes["planet_conf"]
    p_ref = const.reference_surface_pressure
    p_tgt = p_ref_frac * p_ref
    pres = cubelist.extract_cube(model.pres)
    p_tgt.convert_units(pres.units)
    out = interp_to_pres_lev(
        cubelist, constraint, [p_tgt.data], interpolator=interpolator, model=model
    )
    return out


def isel(obj, coord, idx, skip_not_found=None):
    """
    Slice cubes by an index of a coordinate (index-selector).

    Parameters
    ----------
    obj: iris.cube.Cube or iris.cube.CubeList
        Cube or list of cubes.
    coord: str or iris.coords.Coord
        Coordinate for selection.
    idx: int
        Index along the given coordinate.
    skip_not_found: bool or str, optional
        Skip if coordinate not found. By default it is active when dealing with
        a cube list and inactive if dealing with a single cube.

    Returns
    -------
    iris.cube.Cube or iris.cube.CubeList
        Slice along the coordinate.
    """
    if isinstance(obj, iris.cube.Cube):
        try:
            _coord = obj.coord(coord)
            val = _coord.points[idx]
            try:
                val = _coord.units.num2date(val)
            except ValueError:
                pass
            constr = iris.Constraint(**{_coord.name(): lambda x: x.point == val})
            out = obj.extract(constr)
        except CoNotFound as e:
            if skip_not_found:
                out = obj
            else:
                raise e
    elif isinstance(obj, (list, set, tuple)):
        out = []
        for cube in obj:
            out.append(
                isel(cube, coord, idx, skip_not_found=((skip_not_found is None) or skip_not_found))
            )
        out = obj.__class__(out)
    return out


def nearest_coord_value(cube, coord_name, val):
    """
    Get the nearest value of a coordinate.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with the coordinate
    coord_name: str or iris.coords.Coord
        Coordinate where to look the nearest point up
    val: int or float
        The value to find

    Returns
    -------
    int or float
        Element of the coordinate array closest to the given `val`.
    """
    coord = cube.coord(coord_name)
    i = coord.nearest_neighbour_index(val)
    return coord.points[i]


def not_equal_coord_axes(cube1, cube2):
    """Given 2 cubes, return axes of unequal dimensional coordinates."""
    coord_comp = iris.analysis.coord_comparison(cube1, cube2)
    neq_dim_coords = set(coord_comp["not_equal"]).intersection(set(coord_comp["dimensioned"]))
    dims = []
    for coord_pair in neq_dim_coords:
        for coord in coord_pair:
            dims.append(iris.util.guess_coord_axis(coord))
    return set(filter(None, dims))


def regrid_3d(cube, target, model=um):
    """
    Regrid a cube in the horizontal and in the vertical on to coordinates of the target cube.

    Adapted from https://github.com/LSaffin/iris-extensions

    Parameters
    ----------
    cube: iris.cube.Cube
        The cube to be regridded.
    target: iris.cube.Cube
        The cube to regrid to.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.

    Returns
    -------
        iris.cube.Cube
    """
    neq_axes = not_equal_coord_axes(cube, target)
    if neq_axes.intersection(["X", "Y"]):
        cube = cube.regrid(target, iris.analysis.Linear())

    # Interpolate in the vertical if needed
    if "Z" in neq_axes:
        vert_coord = model.z
        if vert_coord is None:
            z = get_dim_coord(target, "Z")
        else:
            z = target.coord(vert_coord)
        cube = cube.interpolate([(z.name(), z.points)], iris.analysis.Linear())
        ensure_bounds(cube, coords=[z])

    # Match coordinate information
    # XXX is this needed?
    # newcube = target.copy(data=cube.data)
    # newcube.rename(cube.name())
    # newcube.units = cube.units

    # Put back correct time information
    # for coord in newcube.aux_coords:
    #     if iris.util.guess_coord_axis(coord) == "T":
    #         newcube.remove_coord(coord)
    # for coord in cube.aux_coords:
    #     if iris.util.guess_coord_axis(coord) == "T":
    #         newcube.add_aux_coord(coord)
    return cube


def replace_z_coord(cube, model=um):
    """
    Replace dimensional vertical coordinate.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.

    Returns
    -------
    iris.cube.Cube
        Copy of the input cube with a new vertical coordinate.
    """
    new_cube = cube.copy()
    new_cube.coord(model.z).bounds = None
    iris.util.promote_aux_coord_to_dim_coord(new_cube, model.z)
    ensure_bounds(new_cube, coords=[model.z])
    for coord in [model.s, model.lev]:
        # By default, model levels and sigma coordinates are removed.
        try:
            new_cube.remove_coord(coord)
        except CoNotFound:
            pass
    return new_cube


def roll_cube_0_360(cube_in, model=um):
    """
    Take a cube spanning -180...180 degrees in longitude and roll it to 0...360 degrees.

    Works with global model output, and in some cases for regional.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    model: aeolus.model.Model, optional
        Model class with a relevant longitude coordinate name.

    Returns
    -------
    iris.cube.Cube

    See also
    --------
    aeolus.coord.roll_cube_pm180
    """
    cube = cube_in.copy()
    coord_name = model.x  # get the name of the longitude coordinate
    lon = cube.coord(coord_name)
    if (lon.points < 0.0).any():
        add = 180
        cube.data = np.roll(cube.data, len(lon.points) // 2, axis=-1)
        if lon.has_bounds():
            bounds = lon.bounds + add
        else:
            bounds = None
        cube.replace_coord(lon.copy(points=lon.points + add, bounds=bounds))
    return cube


def roll_cube_pm180(cube_in, model=um):
    """
    Take a cube spanning 0...360 degrees in longitude and roll it to -180...180 degrees.

    Works with global model output, and in some cases for regional.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    model: aeolus.model.Model, optional
        Model class with a relevant longitude coordinate name.

    Returns
    -------
    iris.cube.Cube

    See also
    --------
    aeolus.coord.roll_cube_0_360
    """
    cube = cube_in.copy()
    coord_name = model.x  # get the name of the longitude coordinate
    xcoord = cube.coord(coord_name)
    if (xcoord.points >= 0.0).all():
        assert is_regular(xcoord), "Operation is only valid for a regularly spaced coordinate."
        if _is_longitude_global(xcoord.points):
            # Shift data symmetrically only when dealing with global cubes
            cube.data = np.roll(cube.data, len(xcoord.points) // 2, axis=-1)

        if xcoord.has_bounds():
            bounds = wrap_lons(xcoord.bounds, -180, 360)  # + subtract
            bounds = bounds[bounds[:, 0].argsort(axis=0)]
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
    return cube


def vertical_cross_section_area(cube2d, r_planet=None):
    """Create a cube of vertical cross-section areas in metres."""
    cube2d = cube2d.copy()
    if r_planet is None:
        r = get_planet_radius(cube2d)
    else:
        r = r_planet
    m_per_deg = (np.pi / 180) * r
    if iris.util.guess_coord_axis(cube2d.dim_coords[1]) == "X":
        m_per_deg *= np.cos(np.deg2rad(cube2d.coord(axis="Y").points[0]))

    for dim_coord in cube2d.dim_coords:
        if not dim_coord.has_bounds():
            dim_coord.guess_bounds()
    x_bounds = cube2d.coord(cube2d.dim_coords[1]).bounds
    z_bounds = cube2d.coord(cube2d.dim_coords[0]).bounds

    vc_area = cube2d.copy(
        data=(z_bounds[:, 1] - z_bounds[:, 0])[:, None]
        * ((x_bounds[:, 1] - x_bounds[:, 0])[None, :] * m_per_deg)
    )
    vc_area.units = "m**2"
    vc_area.rename("vertical_section_area")
    for dim_coord in vc_area.dim_coords:
        dim_coord.bounds = None
    return vc_area


def volume_weights_cube(cube, r_planet=None, normalize=False, model=um):
    """
    Create a cube of volume weights from a grid of a given cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude, latitude and height coordinates
    r_planet: float, optional
        Radius of the planet (m). If not given, an attempt is made
        to get it from the cube metadata.
    normalize: bool, optional
        Normalize the data.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.

    Returns
    -------
    iris.cube.Cube
        Cube of area weights with the same metadata as the input cube
    """
    area_cube = area_weights_cube(cube, r_planet=r_planet, normalize=normalize, model=model)
    height_deltas = coord_delta_to_cube(cube, model.z, normalize=normalize)
    volume = area_cube * height_deltas
    if normalize:
        volume.rename("normalized_volume_weights")
        volume.convert_units("1")
    else:
        volume.rename("volume_weights")
        volume.convert_units("m**3")
    return volume
