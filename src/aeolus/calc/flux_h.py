# -*- coding: utf-8 -*-
"""Integrated fluxes."""
import warnings

import iris

import numpy as np

from ..const import get_planet_radius
from ..coord import nearest_coord_value, vertical_cross_section_area
from ..exceptions import AeolusWarning
from ..model import um


__all__ = ("horizontal_fluxes_through_region_boundaries", "net_horizontal_flux_to_region")


def horizontal_fluxes_through_region_boundaries(
    scalar_cube, region, u, v, r_planet=None, vertical_constraint=None, warn_thresh=10, model=um
):
    """Calculate horizontal fluxes of `scalar_cube` through planes of a rectangular region."""
    perpendicular_wind_cmpnts = {um.x: u, um.y: v}

    if r_planet is None:
        r = get_planet_radius(scalar_cube)
    else:
        r = r_planet

    total_h_fluxes = iris.cube.CubeList()
    for bound in region:
        this_coord = bound["coord"]
        other_coord, (other_min, other_max) = region._perpendicular_side_limits(bound["name"])
        nearest = nearest_coord_value(scalar_cube, this_coord, bound["value"])
        if abs(nearest - bound["value"]) >= warn_thresh:
            warnings.warn(
                f"Nearest value is {np.round(nearest - bound['value'], 2)} deg away"
                f" from the given value of {this_coord}",
                AeolusWarning,
            )
        vcross_cnstr = iris.Constraint(**{this_coord: nearest})
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
        vcross_area = vertical_cross_section_area(cube_slice, r_planet=r)

        # Calculate side flux (2d)
        cube = perpendicular_wind_cmpnts[this_coord].extract(vcross_cnstr) * cube * vcross_area
        cube.rename(f"{scalar_cube.name()}_flux_through_{bound['name']}_boundary")
        # Total flux
        collapsible_dims = [
            i for i in cube.dim_coords if iris.util.guess_coord_axis(i) in ["Z", "Y", "X"]
        ]
        cube_total = cube.collapsed(collapsible_dims, iris.analysis.SUM)
        total_h_fluxes.append(cube_total)
    return total_h_fluxes


def net_horizontal_flux_to_region(
    scalar_cube, region, u, v, r_planet=None, vertical_constraint=None, model=um
):
    """Calculate horizontal fluxes of `scalar_cube` quantity and add them to get the net result."""
    total_h_fluxes = horizontal_fluxes_through_region_boundaries(
        scalar_cube,
        region,
        u,
        v,
        r_planet=r_planet,
        vertical_constraint=vertical_constraint,
        model=model,
    )
    net_flux = (
        total_h_fluxes.extract_cube(
            iris.Constraint(cube_func=lambda x: "through_west" in x.name())
        )
        - total_h_fluxes.extract_cube(
            iris.Constraint(cube_func=lambda x: "through_east" in x.name())
        )
        + total_h_fluxes.extract_cube(
            iris.Constraint(cube_func=lambda x: "through_south" in x.name())
        )
        - total_h_fluxes.extract_cube(
            iris.Constraint(cube_func=lambda x: "through_north" in x.name())
        )
    )
    net_flux.rename(f"net_{scalar_cube.name()}_horizontal_flux_to_region")
    net_flux.attributes["region_str"] = str(region)

    return net_flux
