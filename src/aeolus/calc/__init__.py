"""Science calculations."""
import functools

import iris

from .calculus import d_dx, d_dy, d_dz, deriv, div_h, integrate
from .diag import (
    air_density,
    air_potential_temperature,
    air_temperature,
    bond_albedo,
    dry_lapse_rate,
    flux,
    geopotential_height,
    ghe_norm,
    heat_redist_eff,
    horiz_wind_cmpnts,
    precip_sum,
    sfc_net_energy,
    sfc_water_balance,
    superrotation_index,
    toa_cloud_radiative_effect,
    toa_eff_temp,
    toa_net_energy,
    water_path,
)
from .flux_h import horizontal_fluxes_through_region_boundaries, net_horizontal_flux_to_region
from .stats import (
    cumsum,
    last_n_day_mean,
    meridional_mean,
    minmaxdiff,
    normalize_cube,
    region_mean_diff,
    spatial,
    spatial_mean,
    spatial_quartiles,
    time_mean,
    vertical_mean,
    zonal_mean,
)

__all__ = (
    "air_density",
    "air_potential_temperature",
    "air_temperature",
    "bond_albedo",
    "cumsum",
    "d_dx",
    "d_dy",
    "d_dz",
    "deriv",
    "div_h",
    "dry_lapse_rate",
    "flux",
    "geopotential_height",
    "ghe_norm",
    "heat_redist_eff",
    "horiz_wind_cmpnts",
    "horizontal_fluxes_through_region_boundaries",
    "integrate",
    "last_n_day_mean",
    "meridional_mean",
    "minmaxdiff",
    "net_horizontal_flux_to_region",
    "normalize_cube",
    "precip_sum",
    "region_mean_diff",
    "region_mean_diff",
    "sfc_net_energy",
    "sfc_water_balance",
    "spatial",
    "spatial_mean",
    "spatial_quartiles",
    "superrotation_index",
    "time_mean",
    "toa_cloud_radiative_effect",
    "toa_eff_temp",
    "toa_net_energy",
    "vertical_mean",
    "water_path",
    "zonal_mean",
)


def update_metadata(name=None, units=None):
    """Update metadata of a cube returned by a function."""

    def update_name_units(cube):
        """Update name and convert units."""
        cube.rename(name)
        cube.convert_units(units)

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Call the decorated function
            out = func(*args, **kwargs)
            if isinstance(out, iris.cube.Cube):
                update_name_units(out)
            elif isinstance(out, iris.cube.CubeList):
                [update_name_units(cube) for cube in out]
            return out

        return wrapper

    return decorator
