# -*- coding: utf-8 -*-
"""Model-specific dictionaries of variable names and coordinates."""
from .base import Model

__all__ = ("um", "um_stash")


um = Model(
    # Coordinates
    t="time",
    fcst_ref="forecast_reference_time",
    fcst_prd="forecast_period",
    z="level_height",
    lev="model_level_number",
    s="sigma",
    d="depth",
    p="air_pressure",
    y="latitude",
    x="longitude",
    # Variables
    u="x_wind",
    v="y_wind",
    w="upward_air_velocity",
    pres="air_pressure",
    thta="air_potential_temperature",
    exner="dimensionless_exner_function",
    sh="specific_humidity",
    t_sfc="surface_temperature",
    p_sfc="surface_air_pressure",
    toa_isr="toa_incoming_shortwave_flux",
    toa_olr="toa_outgoing_longwave_flux",
    toa_olr_cs="toa_outgoing_longwave_flux_assuming_clear_sky",
    toa_osr="toa_outgoing_shortwave_flux",
    toa_osr_cs="toa_outgoing_shortwave_flux_assuming_clear_sky",
    sfc_dn_lw="surface_downwelling_longwave_flux_in_air",
    sfc_dn_lw_cs="surface_downwelling_longwave_flux_in_air_assuming_clear_sky",
    sfc_dn_sw="surface_downwelling_shortwave_flux_in_air",
    sfc_dn_sw_cs="surface_downwelling_shortwave_flux_in_air_assuming_clear_sky",
    sfc_net_down_lw="surface_net_downward_longwave_flux",
    sfc_net_down_sw="surface_net_downward_shortwave_flux",
    lw_up="upwelling_longwave_flux_in_air",
    lw_up_forcing="upwelling_longwave_flux_in_air_with_forcing",
    lw_dn_forcing="downwelling_longwave_flux_in_air_with_forcing",
    sw_up="upwelling_shortwave_flux_in_air",
    lw_dn="downwelling_longwave_flux_in_air",
    sw_dn="downwelling_shortwave_flux_in_air",
    lw_up_cs="upwelling_longwave_flux_in_air_assuming_clear_sky",
    sw_up_cs="upwelling_shortwave_flux_in_air_assuming_clear_sky",
    lw_dn_cs="downwelling_longwave_flux_in_air_assuming_clear_sky",
    sw_dn_cs="downwelling_shortwave_flux_in_air_assuming_clear_sky",
    sfc_shf="surface_upward_sensible_heat_flux",
    sfc_lhf="surface_upward_latent_heat_flux",
    sfc_evap="surface_upward_water_flux",
    sfc_x_wind_stress="surface_downward_eastward_stress",
    sfc_y_wind_stress="surface_downward_northward_stress",
    temp="air_temperature",
    dens="air_density",
    ghgt="geopotential_height",
    rh="relative_humidity",
    cld_ice_mf="mass_fraction_of_cloud_ice_in_air",
    cld_liq_mf="mass_fraction_of_cloud_liquid_water_in_air",
    rain_mf="mass_fraction_of_rain_in_air",
    cld_ice_v="ice_cloud_volume_fraction_in_atmosphere_layer",
    cld_liq_v="liquid_cloud_volume_fraction_in_atmosphere_layer",
    cld_v="cloud_volume_fraction_in_atmosphere_layer",
    caf="cloud_area_fraction_assuming_maximum_random_overlap",
    caf_h="high_type_cloud_area_fraction",
    caf_l="low_type_cloud_area_fraction",
    caf_m="medium_type_cloud_area_fraction",
    caf_vl="m01s09i202",
    ppn="precipitation_flux",
    ls_rain="stratiform_rainfall_flux",
    ls_snow="stratiform_snowfall_flux",
    cv_rain="convective_rainfall_flux",
    cv_snow="convective_snowfall_flux",
    dt_sw="change_over_time_in_air_temperature_due_to_shortwave_heating",
    dt_sw_cs="tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky",
    dt_lw="change_over_time_in_air_temperature_due_to_longwave_heating",
    dt_lw_cs="tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky",
    dt_lsppn="change_over_time_in_air_temperature_due_to_stratiform_precipitation",
    dt_bl="change_over_time_in_air_temperature_due_to_boundary_layer_mixing",
    dt_cv="change_over_time_in_air_temperature_due_to_convection",
    dt_lscld="m01s09i181",
    dt_adv="change_over_time_in_air_temperature_due_to_advection",
    dq_sw="change_over_time_in_specific_humidity_due_to_shortwave_heating",
    dq_lw="change_over_time_in_specific_humidity_due_to_longwave_heating",
    dq_lsppn="change_over_time_in_specific_humidity_due_to_stratiform_precipitation",
    dq_bl="change_over_time_in_specific_humidity_due_to_boundary_layer_mixing",
    dq_cv="change_over_time_in_specific_humidity_due_to_convection",
    dq_lscld="m01s09i182",
    dq_adv="change_over_time_in_specific_humidity_due_to_advection",
    du_bl="change_over_time_in_x_wind_due_to_boundary_layer_mixing",
    du_cv="change_over_time_in_x_wind_due_to_convection",
    du_solve="change_over_time_in_x_wind_due_to_pressure_solver",
    du_adv="change_over_time_in_x_wind_due_to_advection",
    du_total="change_over_time_in_x_wind",
    soil_moist="moisture_content_of_soil_layer",
    # Eliassen-Palm flux
    ep_flux_x="northward_eliassen_palm_flux_in_air",
    ep_flux_y="upward_eliassen_palm_flux_in_air",
    ep_flux_div="tendency_of_eastward_wind_due_to_eliassen_palm_flux_divergence",
)

um_stash = Model(
    # Coordinates
    t="time",
    fcst_ref="forecast_reference_time",
    fcst_prd="forecast_period",
    z="level_height",
    lev="model_level_number",
    s="sigma",
    d="depth",
    p="air_pressure",
    y="latitude",
    x="longitude",
    # Variables
    u="m01s00i002",
    v="m01s00i003",
    w="m01s00i150",
    pres="m01s00i408",
    thta="m01s00i004",
    exner="m01s00i255",
    sh="m01s00i010",
    t_sfc="m01s00i024",
    p_sfc="m01s00i409",
    toa_isr="m01s01i207",
    toa_olr="m01s02i205",
    toa_olr_cs="m01s02i206",
    toa_osr="m01s01i208",
    toa_osr_cs="m01s01i209",
    sfc_dn_lw="m01s02i207",
    sfc_dn_lw_cs="m01s02i208",
    sfc_dn_sw="m01s01i235",
    sfc_dn_sw_cs="m01s01i210",
    sfc_net_down_lw="m01s02i201",
    sfc_net_down_sw="m01s01i201",
    lw_up="m01s02i217",
    lw_up_forcing="m01s02i417",
    lw_dn_forcing="m01s02i418",
    sw_up="m01s01i217",
    lw_dn="m01s02i218",
    sw_dn="m01s01i218",
    lw_up_cs="m01s02i219",
    sw_up_cs="m01s01i219",
    lw_dn_cs="m01s02i220",
    sw_dn_cs="m01s01i220",
    sfc_shf="m01s03i217",
    sfc_lhf="m01s03i234",
    temp="m01s03i236",
    rh="m01s30s113",
    cld_ice_mf="m01s00i012",
    cld_liq_mf="m01s00i254",
    rain_mf="m01s00i272",
    cld_ice_v="m01s00i268",
    cld_liq_v="m01s00i267",
    cld_v="m01s00i266",
    caf="m01s09i217",
    caf_h="m01s09i205",
    caf_m="m01s09i204",
    caf_l="m01s09i203",
    caf_vl="m01s09i202",
    ppn="m01s05i216",
    ls_rain="m01s04i203",
    ls_snow="m01s04i204",
    cv_rain="m01s05i205",
    cv_snow="m01s05i206",
    dt_sw="m01s01i181",
    dt_sw_cs="m01s01i233",
    dt_lw="m01s02i181",
    dt_lw_cs="m01s02i233",
    dt_lsppn="m01s04i181",
    dt_bl="m01s03i181",
    dt_cv="m01s05i181",
    dt_lscld="m01s09i181",
    dt_adv="m01s12i181",
    dq_sw="m01s01i182",
    dq_lw="m01s02i182",
    dq_lsppn="m01s04i182",
    dq_bl="m01s03i182",
    dq_cv="m01s05i182",
    dq_lscld="m01s09i182",
    dq_adv="m01s12i182",
    soil_moist="m01s08i223",
    # Eliassen-Palm flux
    ep_flux_x="m01s30i312",
    ep_flux_y="m01s30i313",
    ep_flux_div="m01s30i314",
    du_bl="m01s03i185",
    du_cv="m01s05i185",
    du_solve="m01s10i185",
    du_adv="m01s12i185",
    du_total="m01s30i185",
)
