"""Utilities for working with LFRic output."""

from collections.abc import Sequence
from pathlib import Path
from typing import Callable, Optional, Union

import iris
import iris.coords
from iris.cube import Cube, CubeList
import iris.exceptions
from iris.mesh import MeshCoord, MeshXY
import iris.util
import numpy as np

from .io import load_vert_lev
from .model import lfric
from .model.base import Model

__all__ = (
    "add_um_height_coord",
    "add_equally_spaced_height_coord",
    "clean_attrs",
    "fix_time_coord",
    "load_lfric_raw",
    "replace_level_coord_with_height",
    "replace_mesh",
    "simple_regrid_lfric",
    "ugrid_spatial",
    "ugrid_spatial_mean",
)


def add_um_height_coord(
    cube: Cube, field: str, filename: str, path_to_levels_file: Path
) -> Cube:
    """
    Callback for `iris.load` specifically for LFRic data with vertical coord.

    Adds a vertical coordinate of level heights from a UM vertical levels file.

    Parameters
    ----------
    path_to_levels_file: Path-like
        Path to the vertical levels file.

    Examples
    --------
    >>> from functools import partial
    >>> iris.load(
        filename,
        callback=partial(
            add_um_height_coord,
            path_to_levels_file=Path('/path/to/vertlevs_L38_29t_9s_40km')
        )
    )
    """
    lev_coords = [
        i.name() for i in cube.dim_coords if i.name().endswith("_levels")
    ]
    if len(lev_coords) == 1:
        thlev = load_vert_lev(path_to_levels_file, lev_type="theta")
        rholev = load_vert_lev(path_to_levels_file, lev_type="rho")
        lev_coord = cube.coord(lev_coords[0])
        if "full" in lev_coord.name().lower():
            points = thlev
        else:
            points = rholev
        hgt_coord = iris.coords.AuxCoord(
            points,
            long_name="level_height",
            units="m",
            attributes={"positive": "up"},
        )
        hgt_coord.guess_bounds()
        cube.add_aux_coord(hgt_coord, data_dims=cube.coord_dims(lev_coord))

    return cube


def add_equally_spaced_height_coord(
    cube: Cube, field: str, filename: str, model_top_height: float
) -> Cube:
    """
    Callback for `iris.load` specifically for LFRic data with vertical coord.

    Adds a vertical coordinate of level heights assuming equally spaced levels.

    Examples
    --------
    >>> from functools import partial
    >>> iris.load(
        filename,
        callback=partial(
            add_equally_spaced_height_coord,
            model_top_height=32000
        )
    )
    """
    lev_coords = [
        i.name() for i in cube.dim_coords if i.name().endswith("_levels")
    ]
    if len(lev_coords) == 1:
        lev_coord = cube.coord(lev_coords[0])
        if "full" in lev_coord.name().lower():
            points = np.linspace(0, model_top_height, lev_coord.shape[0])
        else:
            points = np.linspace(0, model_top_height, lev_coord.shape[0] + 1)
            points = 0.5 * (points[:-1] + points[1:])
        hgt_coord = iris.coords.AuxCoord(
            points,
            long_name="level_height",
            units="m",
            attributes={"positive": "up"},
        )
        hgt_coord.guess_bounds()
        cube.add_aux_coord(hgt_coord, data_dims=cube.coord_dims(lev_coord))

    return cube


def clean_attrs(cube: Cube, field: str, filename: str) -> Cube:
    """
    Callback for `iris.load` to clean up trivial attributes in LFRic data.

    Needed for concatenating cubes.
    """
    for attr in ["timeStamp", "uuid"]:
        try:
            cube.attributes.pop(attr)
        except KeyError:
            pass
    return cube


def fix_time_coord(
    cube: Cube,
    field: str,
    filename: str,
    downgrade_to_scalar: Optional[bool] = False,
) -> Cube:
    """
    Callback function for `iris.load` specifically for UGRID data.

    All field variables should have time as a coord, and mesh variables
    will not have a time coord. So ignore anything that doesn't have a
    time coord. Note that mesh info will still be loaded in the
    coordinate part of the cube.

    Make sure that time coordinate is correct type, seems to always be
    loaded as a basic AuxCoord rather than a Scalar or DimCoord.
    """
    # Ignore cubes with no time coordinate
    if len(cube.coords("time")) == 0:
        raise iris.exceptions.IgnoreCubeException

    # Else, fix metadata in variables
    else:
        tcoord = cube.coord("time")
        # If only 1 time coordinate value, downgrade AuxCoord to Scalar
        if downgrade_to_scalar and tcoord.points.size == 1:
            cube = cube.extract(iris.Constraint(time=tcoord.cell(0)))
        # Else, upgrade AuxCoord to DimCoord
        else:
            iris.util.promote_aux_coord_to_dim_coord(cube, tcoord)
    return cube


def load_lfric_raw(
    fnames: Sequence[Union[Path, str]],
    callback: Optional[Callable] = None,
    drop_coord: Optional[Sequence[str]] = (),
) -> CubeList:
    """Load raw LFRic data."""
    cl_raw = iris.load(fnames, callback=callback)
    cl_raw = CubeList(clean_attrs(cube, None, None) for cube in cl_raw)
    for coord in drop_coord:
        for cube in cl_raw:
            try:
                cube.remove_coord(coord)
            except iris.exceptions.CoordinateNotFoundError:
                pass
    cl_raw = cl_raw.concatenate(check_aux_coords=False)
    return cl_raw


def replace_level_coord_with_height(cube: Cube) -> Cube:
    """
    Remove full_ or half_levels coordinate and replace it with level_height.

    Assumes that the cube contains level_height already.

    See Also
    --------
    add_equally_spaced_height_coord
    """
    lev_coords = [
        i.name() for i in cube.dim_coords if i.name().endswith("_levels")
    ]
    if len(lev_coords) == 1:
        lev_coord = cube.coord(lev_coords[0])
        cube.remove_coord(lev_coord)
        iris.util.promote_aux_coord_to_dim_coord(cube, "level_height")
    return cube


def replace_mesh(cube: Cube, new_mesh: MeshXY) -> Cube:
    """Replace mesh in a 1d cube by creating a new copy of that cube."""
    mesh_x = MeshCoord(mesh=new_mesh, location="face", axis="x")
    mesh_y = MeshCoord(mesh=new_mesh, location="face", axis="y")

    new_cube = Cube(cube.data, aux_coords_and_dims=[(mesh_x, 0), (mesh_y, 0)])
    new_cube.metadata = cube.metadata
    return new_cube


def reshape_ugrid_cube_to_xy_grid(
    cube: iris.cube.Cube, nx: int, ny: int, delta_x: float, delta_y: float
) -> iris.cube.Cube:
    """
    Reshape an iris cube on a UGRID into a cube with x/y distance coordinates.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube to reshape
    nx : int, optional
        Number of grid points in x direction
    ny : int, optional
        Number of grid points in y direction
    delta_x : float, optional
        Grid spacing in x direction in metres
    delta_y : float, optional
        Grid spacing in y direction in metres

    Returns
    -------
    iris.cube.Cube
        Reshaped cube with x/y coordinates
    """
    # Promote auxiliary time coordinate to dimension coordinate if necessary
    try:
        iris.util.promote_aux_coord_to_dim_coord(cube, "time")
    except iris.exceptions.CoordinateNotFoundError:
        pass

    # Calculate domain sizes
    domain_size_x: float = nx * delta_x
    domain_size_y: float = ny * delta_y

    # Define new shape
    new_shape: tuple = cube.shape[:-1] + (ny, nx)

    # Create spatial coordinates
    x_coord: iris.coords.DimCoord = iris.coords.DimCoord(
        np.linspace(0, domain_size_x, nx + 1)[:-1],
        standard_name="projection_x_coordinate",
        units="m",
    )
    y_coord: iris.coords.DimCoord = iris.coords.DimCoord(
        np.linspace(0, domain_size_y, ny + 1)[:-1],
        standard_name="projection_y_coordinate",
        units="m",
    )

    # Create new cube with reshaped data
    new_cube: iris.cube.Cube = iris.cube.Cube(
        iris.util.reverse(
            cube.core_data().reshape(new_shape),
            new_shape.index(y_coord.shape[0]),
        ),
        dim_coords_and_dims=[
            (c, cube.coord_dims(c)[0]) for c in cube.dim_coords
        ]
        + [
            (y_coord, len(new_shape) - 2),
            (x_coord, len(new_shape) - 1),
        ],
    )

    # Copy metadata from original cube
    new_cube.metadata = iris.common.metadata.CubeMetadata(
        **cube.metadata._asdict()
    )

    return new_cube


def _regrid_cubes_on_edges(
    src_cube: Cube, tgt_cube: Cube, method: Optional[str] = "auto"
) -> Cube:
    """Temporary solution to regrid cubes defined on cell edges."""
    # Using the approach in gungho/rose-stem/app/plot/bin/plot_zonal.py
    from iris._concatenate import _CubeSignature
    from iris.analysis.cartography import get_xy_grids
    from scipy.interpolate import griddata

    src_x = src_cube.coord("longitude").points
    src_y = src_cube.coord("latitude").points
    tgt_x2d, tgt_y2d = get_xy_grids(tgt_cube)

    non_xy_coords = [
        name
        for c in src_cube.dim_coords
        if (name := c.name()) not in ["longitude", "latitude"]
    ]

    result = CubeList()
    for src_2d_slice in src_cube.slices_over(non_xy_coords):
        # Interpolate using delaunay triangularization
        if method == "auto":
            reg_data_linear = griddata(
                (src_x, src_y),
                src_2d_slice.data,
                (tgt_x2d, tgt_y2d),
                method="linear",
            )
            reg_data_nearest = griddata(
                (src_x, src_y),
                src_2d_slice.data,
                (tgt_x2d, tgt_y2d),
                method="nearest",
            )
            reg_data_linear[np.isnan(reg_data_linear)] = reg_data_nearest[
                np.isnan(reg_data_linear)
            ]
            reg_data = reg_data_linear
        else:
            reg_data = griddata(
                (src_x, src_y),
                src_2d_slice.data,
                (tgt_x2d, tgt_y2d),
                method=method,
            )
        reg_cube_slice = tgt_cube.copy(data=reg_data)
        for scalar_coord in _CubeSignature(src_2d_slice).scalar_coords:
            reg_cube_slice.add_aux_coord(scalar_coord)
        result.append(reg_cube_slice)
    result = result.merge_cube()
    dim_order = [*range(src_cube.ndim + 1)]
    dim_order[0], dim_order[1] = 1, 0
    result.transpose(dim_order)
    result.metadata = src_cube.metadata
    return result


def simple_regrid_lfric(
    cube_list: CubeList,
    tgt_cube: Cube,
    ref_cube_constr: Optional[str] = "air_potential_temperature",
    interp_vertically: Optional[bool] = True,
    model: Optional[Model] = lfric,
) -> CubeList:
    """Quick&dirty regridding of LFRic data to a common height/lat/lon grid."""
    from esmf_regrid.experimental.unstructured_scheme import (
        MeshToGridESMFRegridder,
    )

    # Horizontal regridding
    result_h = CubeList()
    ref_cube = cube_list.extract_cube(ref_cube_constr)
    regridder = MeshToGridESMFRegridder(
        ref_cube, tgt_cube, method="conservative"
    )
    for cube in cube_list:
        if cube.location == "edge":
            cube_reg = _regrid_cubes_on_edges(cube, tgt_cube)
        else:
            cube_reg = regridder(cube)
        result_h.append(cube_reg)
    result = result_h
    if interp_vertically:
        # Vertical interpolation
        result_v = iris.cube.CubeList()
        tgt_points = (model.z, ref_cube.coord(model.z).points)
        for cube in result:
            if len(
                [
                    i.name()
                    for i in cube.dim_coords
                    if i.name().endswith("_levels")
                ]
            ):
                result_v.append(
                    replace_level_coord_with_height(cube).interpolate(
                        [tgt_points], iris.analysis.Linear()
                    )
                )
            else:
                result_v.append(cube)
        result = result_v
    return result


def ugrid_spatial(
    cube: Cube,
    aggr: str,
    model: Optional[Model] = lfric,
    **kwargs,
) -> Cube:
    """Collapse a UGRID `iris.cube.Cube` over x and y coord with no weights."""
    cube_copy = cube.copy()
    tmp_coord = iris.coords.AuxCoord(
        points=np.arange(cube.coord(model.x).shape[0]),
        long_name="mesh_coordinates",
    )

    cube_copy.add_aux_coord(tmp_coord, data_dims=cube.coord_dims(model.x))

    cube_copy.remove_coord(model.x)
    cube_copy.remove_coord(model.y)

    cube_aggr = cube_copy.collapsed(
        tmp_coord.name(), getattr(iris.analysis, aggr.upper()), **kwargs
    )
    cube_aggr.remove_coord(tmp_coord)
    return cube_aggr


def ugrid_spatial_mean(
    cube: Cube, model: Optional[Model] = lfric, **kwargs
) -> Cube:
    """Average a UGRID `iris.cube.Cube` spatially."""
    return ugrid_spatial(cube, "mean", model=model, **kwargs)
