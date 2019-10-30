"""Science calculations."""
from .calculus import integrate
from .diag import sfc_water_balance, toa_cloud_radiative_effect, toa_net, total_precip
from .stats import last_year_mean, meridional_mean, spatial, spatial_quartiles

__all__ = (
    "integrate",
    "sfc_water_balance",
    "toa_cloud_radiative_effect",
    "toa_net",
    "total_precip",
    "last_year_mean",
    "meridional_mean",
    "spatial",
    "spatial_quartiles",
)
