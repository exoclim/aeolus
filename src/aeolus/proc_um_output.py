# -*- coding: utf-8 -*-
"""Functions for post-processing UM output."""
import warnings
from typing import Optional

import iris
from iris.cube import CubeList

from .coord import (
    add_planet_calendar,
    regrid_3d,
    replace_z_coord,
    roll_cube_pm180,
    ensure_bounds,
    interp_to_cube_time,
)
from .model import um
from .model.base import Model
from .subset import CM_INST_CONSTR, CM_MEAN_CONSTR, DimConstr, unique_cubes


__all__ = ("process_cubes",)


def process_cubes(
    cubelist: CubeList,
    timestep: Optional[int, float] = 1,
    ref_cube_constr: Optional[str] = um.thta,
    extract_incr: Optional[bool] = True,
    extract_mean: Optional[bool] = True,
    regrid_multi_lev: Optional[bool] = True,
    roll_pm180: Optional[bool] = True,
    add_calendar: Optional[bool] = False,
    calendar: Optional[dict] = {},
    planet: Optional[str] = "earth",
    remove_duplicates: Optional[bool] = True,
    use_varpack: Optional[bool] = False,
    varpack: Optional[dict] = {},
    interp_time: Optional[bool] = True,
    model: Optional[Model] = um,
) -> CubeList:
    """Post-process data for easier analysis."""
    DC = DimConstr(model=model)

    if remove_duplicates:
        cubelist = unique_cubes(cubelist)

    cm_constr = iris.Constraint()
    if extract_mean:
        cm_constr &= CM_MEAN_CONSTR
    else:
        cm_constr &= CM_INST_CONSTR
    cubelist = cubelist.extract(cm_constr)

    cubes = CubeList()

    # First, extract all multi-level fields
    if use_varpack:
        cubes += cubelist.extract(varpack["multi_level"])

        # Increments
        if extract_incr:
            cubes += cubelist.extract(
                iris.AttributeConstraint(
                    STASH=lambda x: x.item in [181, 182, 233] and x.section not in [0]
                )
            )
    else:
        cubes = cubelist.extract(DC.strict.tmyx)

    if regrid_multi_lev:
        # Interpolation & regridding to common grid
        ref_cube = cubes.extract_cube(ref_cube_constr)
        ref_cube = replace_z_coord(ref_cube, model=model)

        # Interpolate to common levels
        cubes = CubeList(
            [regrid_3d(replace_z_coord(cube, model=model), ref_cube, model=model) for cube in cubes]
        )
    # Fix units of increments
    for cube in cubes:
        if cube.attributes["STASH"].item in [181, 233]:
            incr_unit = "K"
        elif cube.attributes["STASH"].item in [182]:
            incr_unit = "kg kg^-1"
        if cube.attributes["STASH"].item == 233:
            cube.units = f"{incr_unit} s^-1"
        elif cube.attributes["STASH"].item in [181, 182]:
            cube.units = f"{1/timestep} {incr_unit} s^-1"
            cube.convert_units(f"{incr_unit} s^-1")

    # Add all single-level cubes
    if use_varpack:
        cubes += cubelist.extract(varpack["single_level"])
    else:
        cubes += cubelist.extract(DC.strict.tyx)

    # Roll cubes to +/- 180 degrees in longitude for easier analysis
    if roll_pm180:
        rolled_cubes = CubeList()
        for cube in cubes:
            r_c = roll_cube_pm180(cube)
            ensure_bounds(r_c)
            rolled_cubes.append(r_c)
    else:
        rolled_cubes = cubes

    if interp_time:
        ref_cube = rolled_cubes.extract_cube(model.t_sfc)  # TODO: make an arg
        final = CubeList()
        for cube in rolled_cubes:
            final.append(interp_to_cube_time(cube, ref_cube, model=model))
    else:
        final = rolled_cubes

    if add_calendar:
        try:
            cal = calendar[planet]
            for cube in final:
                add_planet_calendar(
                    cube,
                    model.t,
                    days_in_year=cal["year"],
                    days_in_month=cal["month"],
                    days_in_day=cal["day"],
                    planet=planet,
                )
        except KeyError:
            warnings.warn(f"Calendar for {planet=} not found.")
    return final
