# -*- coding: utf-8 -*-
"""LFRic-specific dictionaries of variable names and coordinates."""
from .base import Model

__all__ = ("lfric",)


lfric = Model(
    # Coordinates
    t="time",
    z="level_height",
    s="atmosphere_hybrid_sigma_pressure_coordinate",
    lev="full_levels",
    y="latitude",
    x="longitude",
    # Variables
    u="u_in_w3",
    v="v_in_w3",
    w="w_in_wth",
    pres="pressure_in_wth",
    thta="theta",
    # exner="exner_pressure",
    exner="exner_in_wth",
    temp="temperature",
    t_sfc="grid_surface_temperature",
    dens="density",
    p_sfc="pmsl",
    ghgt="geopotential_height",
    toa_isr="sw_direct_toa",
    toa_osr="sw_up_toa",
    toa_olr="lw_up_toa",
    toa_osr_cs="sw_up_clear_toa_rts",
    toa_olr_cs="lw_up_clear_toa_rts",
    caf="cloud_amount_maxrnd",
    sfc_dn_sw="sw_down_surf",
    sfc_dn_lw="lw_down_surf",
    sfc_up_lw="lw_up_surf",
    sfc_up_sw="sw_up_surf",
)
