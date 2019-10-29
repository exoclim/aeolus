# -*- coding: utf-8 -*-
"""Miscellaneous."""
import warnings

import iris

import numpy as np

from .coord_utils import nearest_coord_value
from .exceptions import AeolusWarning


__all__ = (
    "vertical_cross_section_area",
    "horizontal_fluxes_through_region_boundaries",
    "net_horizontal_flux_to_region",
)


def vertical_cross_section_area(cube2d, r_planet):
    """Create a cube of vertical cross-section areas in metres."""
    cube2d = cube2d.copy()
    m_per_deg = (np.pi / 180) * r_planet
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


def horizontal_fluxes_through_region_boundaries(
    scalar_cube, region, u, v, r_planet, vertical_constraint=None
):
    """Calculate horizontal fluxes of `scalar_cube` through planes of a rectangular region."""
    perpendicular_wind_cmpnts = {"longitude": u, "latitude": v}

    total_h_fluxes = iris.cube.CubeList()
    for bound in region:
        this_coord = bound["coord"]
        other_coord, (other_min, other_max) = region._perpendicular_side_limits(bound["name"])
        nearest = nearest_coord_value(scalar_cube, this_coord, bound["value"])
        if abs(nearest - bound["value"]) > 10:
            warnings.warn(
                f"Nearest value is {np.round(nearest - bound['value'], 2)} deg away"
                f" from the given value of {this_coord}",
                AeolusWarning,
            )
        vcross_cnstr = iris.Constraint(**{this_coord: nearest})
        if vertical_constraint is not None:
            vcross_cnstr &= vertical_constraint
        if other_max >= other_min:
            vcross_cnstr &= iris.Constraint(**{other_coord: lambda x: other_min <= x <= other_max})
            cube = scalar_cube.extract(vcross_cnstr)
        else:
            vcross_cnstr &= iris.Constraint(
                **{other_coord: lambda x: (other_max >= x) or (other_min <= x)}
            )
            cube = scalar_cube.extract(vcross_cnstr)
        cube_slice = next(cube.slices([cube.coord(axis="z").name(), other_coord]))
        vcross_area = vertical_cross_section_area(cube_slice, r_planet=r_planet)

        # Calculate energy flux (2d)
        cube = (
            perpendicular_wind_cmpnts[this_coord].extract(vcross_cnstr)
            * scalar_cube.extract(vcross_cnstr)
            * vcross_area
        )
        cube.rename(f"{scalar_cube.name()}_flux_through_{bound['name']}_boundary")
        # Total flux
        cube_total = cube.collapsed(cube.dim_coords, iris.analysis.SUM)
        total_h_fluxes.append(cube_total)
    return total_h_fluxes


def net_horizontal_flux_to_region(scalar_cube, region, u, v, r_planet, vertical_constraint=None):
    """Calculate horizontal fluxes of `scalar_cube` quantity and add them to get the net result."""
    total_h_fluxes = horizontal_fluxes_through_region_boundaries(
        scalar_cube, region, u, v, r_planet, vertical_constraint=vertical_constraint
    )
    net_flux = (
        total_h_fluxes.extract_strict(
            iris.Constraint(cube_func=lambda x: "through_west" in x.name())
        )
        - total_h_fluxes.extract_strict(
            iris.Constraint(cube_func=lambda x: "through_east" in x.name())
        )
        + total_h_fluxes.extract_strict(
            iris.Constraint(cube_func=lambda x: "through_south" in x.name())
        )
        - total_h_fluxes.extract_strict(
            iris.Constraint(cube_func=lambda x: "through_north" in x.name())
        )
    )
    net_flux.rename(f"net_{scalar_cube.name()}_horizontal_flux_to_region")
    net_flux.attributes["region_str"] = str(region)

    return net_flux
