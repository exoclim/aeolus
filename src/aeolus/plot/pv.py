"""Prepare data in spherical coordinates for pyvista plotting."""

from typing import Optional

import geovista as gv
import iris
import numpy as np
import pyvista as pv
from pyvista import grid_from_sph_coords, transform_vectors_sph_to_cart

from ..coord import _cell_bounds
from ..exceptions import _warn
from ..model import um
from ..model.base import Model

__all__ = (
    "grid_for_scalar_cube_sph",
    "grid_for_vector_cubes_sph",
    "ugrid_mesh_to_gv_mesh",
)


def grid_for_scalar_cube_sph(
    cube: iris.cube.Cube,
    z_scale: Optional[float] = 1,
    z_offset: Optional[float] = 0,
    grid: Optional[pv.StructuredGrid] = None,
    label: Optional[str] = "scalar3d",
    model: Optional[Model] = um,
) -> pv.PolyData:
    """
    Create a `pyvista` grid for an `iris` cube in spherical coordinates.

    Parameters
    ----------
    cube: iris.cube.Cube
        2D or 3D cube with (longitude, latitude, [level_height]) coordinates.
    z_scale: float, optional
        Scaling factor for the vertical level dimension.
    z_offset: float, optional
        Scaling offset for the vertical level dimension.
    grid: pyvista.StructuredGrid, optional
        If given, add data to the existing grid, otherwise create a new one
        from the input cube's coordinates.
    label: str, optional
        Label for the data within the grid.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    pyvista.StructuredGrid
       PyVista grid with data in the cell_data attribute.
    """
    if grid is None:
        lons = _cell_bounds(cube.coord(model.x).points)
        lats = 90.0 - _cell_bounds(cube.coord(model.y).points)
        if cube.ndim == 3:
            # TODO check if _cell_bounds should be here
            levels = _cell_bounds(cube.coord(model.z).points)
        else:
            levels = np.array([0])
        levels = z_scale * levels + z_offset

        grid = grid_from_sph_coords(lons, lats, levels)

    # Add data arras to grid
    if label in grid.cell_data:
        _warn(f"Label '{label}' exists in {grid}")
    grid.cell_data[label] = np.array(cube.data).swapaxes(-2, -1).ravel("C")
    return grid


def grid_for_vector_cubes_sph(
    u: iris.cube.Cube,
    v: iris.cube.Cube,
    w: iris.cube.Cube,
    vector_scale: Optional[float] = 1,
    vertical_wind_scale: Optional[float] = 1,
    z_scale: Optional[float] = 1,
    z_offset: Optional[float] = 0,
    xstride: Optional[int] = 1,
    ystride: Optional[int] = 1,
    grid: Optional[pv.StructuredGrid] = None,
    label: Optional[str] = "vector3d",
    model: Optional[Model] = um,
) -> pv.PolyData:
    """
    Take winds in spherical coordinates and create a `pyvista` grid for them.

    Parameters
    ----------
    u: iris.cube.Cube
        2D or 3D cube of x-wind component (zonal wind).
    v: iris.cube.Cube
        2D or 3D cube of y-wind component (meridional wind).
    w: iris.cube.Cube
        2D or 3D cube of z-wind component (vertical wind).
    vector_scale: float, optional
        Scaling factor for vectors.
    vertical_wind_scale: float, optional
        Scaling factor for the vertical wind component (for better visibility).
    z_scale: float, optional
        Scaling factor for the vertical level dimension.
    z_offset: float, optional
        Scaling offset for the vertical level dimension.
    xstride: float, optional
        Stride along the longitude axis.
    ystride: float, optional
        Stride along the latitude axis.
    grid: pyvista.StructuredGrid, optional
        If given, add data to the existing grid, otherwise create a new one
        from the input cube's coordinates.
    label: str, optional
        Label for the data within the grid.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    pyvista.StructuredGrid
       PyVista grid with vector data in the point_data attribute.
    """
    assert (u.shape == v.shape) and (
        u.shape == w.shape
    ), "Wind components should have the same array size!"

    lons = u.coord(model.x).points
    lats = 90.0 - u.coord(model.y).points
    levels = z_scale * u.coord(model.z).points + z_offset

    # Stride vectors
    lons_s = lons[::xstride]
    lats_s = lats[::ystride]

    inv_axes = [*range(u.ndim - 1, -1, -1)]

    u_sdata = u.data[..., ::ystride, ::xstride].transpose(inv_axes)
    v_sdata = v.data[..., ::ystride, ::xstride].transpose(inv_axes)
    w_sdata = w.data[..., ::ystride, ::xstride].transpose(inv_axes)

    # Rescale vertical wind for better viz
    w_sdata *= vertical_wind_scale

    # Create a stack of vectors transformed from spherical to cartesian coords
    vectors = np.stack(
        [
            i.transpose(inv_axes).swapaxes(u.ndim - 2, u.ndim - 1).ravel("C")
            for i in transform_vectors_sph_to_cart(
                lons_s, lats_s, levels, u_sdata, -v_sdata, w_sdata
            )
        ],
        axis=1,
    )
    # Scale vectors
    vectors *= vector_scale

    if grid is None:
        grid = grid_from_sph_coords(lons_s, lats_s, levels)

    # Add vectors to the grid
    grid.point_data[label] = vectors
    return grid


def ugrid_mesh_to_gv_mesh(
    topology: iris.experimental.ugrid.mesh.Mesh, **kwargs
) -> pv.PolyData:
    """Create a geovista mesh from the LFRic mesh."""
    face_node = topology.face_node_connectivity

    indices = face_node.indices_by_location()

    lons, lats = topology.node_coords

    return gv.Transform.from_unstructured(
        lons.points,
        lats.points,
        indices,
        start_index=face_node.start_index,
        **kwargs,
    )
