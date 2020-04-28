"""Science calculations."""
from .calculus import integrate
from .diag import (
    bond_albedo,
    ghe_norm,
    heat_redist_eff,
    precip_sum,
    sfc_net_energy,
    sfc_water_balance,
    toa_cloud_radiative_effect,
    toa_eff_temp,
    toa_net_energy,
    water_path,
)
from .flux_h import horizontal_fluxes_through_region_boundaries, net_horizontal_flux_to_region
from .stats import (
    last_year_mean,
    meridional_mean,
    minmaxdiff,
    region_mean_diff,
    spatial,
    spatial_quartiles,
    vertical_cumsum,
    vertical_mean,
    zonal_mean,
)

__all__ = (
    "integrate",
    "bond_albedo",
    "ghe_norm",
    "heat_redist_eff",
    "horizontal_fluxes_through_region_boundaries",
    "minmaxdiff",
    "net_horizontal_flux_to_region",
    "precip_sum",
    "sfc_net_energy",
    "sfc_water_balance",
    "region_mean_diff",
    "toa_cloud_radiative_effect",
    "toa_eff_temp",
    "toa_net_energy",
    "last_year_mean",
    "meridional_mean",
    "zonal_mean",
    "spatial",
    "spatial_quartiles",
    "vertical_cumsum",
    "vertical_mean",
    "water_path",
)
