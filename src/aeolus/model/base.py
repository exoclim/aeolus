"""Model-specific dictionaries of variable names and coordinates."""
from dataclasses import dataclass


__all__ = "Model"


@dataclass
class Model:
    """Base class for model-specific names."""

    # Coordinates
    t: str = None  # time
    fcst_ref: str = None  # forecast reference
    fcst_prd: str = None  # forecast period
    z: str = None  # height
    lev: str = None  # level number
    s: str = None  # sigma
    y: str = None  # latitude
    x: str = None  # longitude
    # Variables
    # Main
    u: str = None
    v: str = None
    w: str = None
    pres: str = None
    thta: str = None
    exner: str = None
    sh: str = None
    t_sfc: str = None
    p_sfc: str = None
    # Radiation
    toa_isr: str = None
    toa_olr: str = None
    toa_olr_cs: str = None
    toa_osr: str = None
    toa_osr_cs: str = None
    sfc_dn_lw: str = None
    sfc_dn_lw_cs: str = None
    sfc_dn_sw: str = None
    sfc_dn_sw_cs: str = None
    sfc_net_down_lw: str = None
    sfc_net_down_sw: str = None
    lw_up: str = None
    sw_up: str = None
    lw_dn: str = None
    sw_dn: str = None
    lw_up_cs: str = None
    sw_up_cs: str = None
    lw_dn_cs: str = None
    sw_dn_cs: str = None
    # BL
    sfc_shf: str = None
    sfc_lhf: str = None
    sfc_evap: str = None
    # Extra physics
    temp: str = None
    dens: str = None
    ghgt: str = None
    rh: str = None
    # Precip & Cloud
    cld_ice_mf: str = None
    cld_liq_mf: str = None
    rain_mf: str = None
    cld_ice_v: str = None
    cld_liq_v: str = None
    cld_v: str = None
    caf: str = None
    caf_h: str = None
    caf_m: str = None
    caf_l: str = None
    caf_vl: str = None
    ppn: str = None
    ls_rain: str = None
    ls_snow: str = None
    cv_rain: str = None
    cv_snow: str = None
    # Increments of temperature
    dt_sw: str = None
    dt_sw_cs: str = None
    dt_lw: str = None
    dt_lw_cs: str = None
    dt_lsppn: str = None
    dt_bl: str = None
    dt_cv: str = None
    dt_lscld: str = None
    dt_adv: str = None
    # Increments of specific humidity
    dq_sw: str = None
    dq_lw: str = None
    dq_lsppn: str = None
    dq_bl: str = None
    dq_cv: str = None
    dq_lscld: str = None
    dq_adv: str = None

    def __repr__(self):
        """Override the repr method of the dataclass."""
        size = len([_ for _, v in self.__dataclass_fields__.items() if v.name is not None])
        return f"{self.__class__}({size} fields)"
