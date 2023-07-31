# -*- coding: utf-8 -*-
"""Some commonly used diagnostics in atmospheric science."""
from cf_units import Unit

from iris.analysis.calculus import _coord_cos
from iris.analysis.maths import add, apply_ufunc, multiply
from iris.exceptions import ConstraintMismatchError as ConMisErr
from iris.util import reverse

import numpy as np

from .calculus import d_dz, integrate
from .stats import cumsum, spatial_mean, time_mean, zonal_mean
from ..const import init_const
from ..coord import coord_to_cube, ensure_bounds, regrid_3d
from ..exceptions import ArgumentError, MissingCubeError
from ..meta import const_from_attrs, preserve_shape, update_metadata
from ..model import um
from ..subset import _dim_constr


__all__ = (
    "air_density",
    "air_potential_temperature",
    "air_pressure",
    "air_temperature",
    "bond_albedo",
    "calc_derived_cubes",
    "dry_lapse_rate",
    "flux",
    "geopotential_height",
    "ghe_norm",
    "greenhouse_effect",
    "heat_redist_eff",
    "horiz_wind_cmpnts",
    "meridional_mass_streamfunction",
    "precip_sum",
    "sfc_net_energy",
    "sfc_water_balance",
    "sigma_p",
    "toa_cloud_radiative_effect",
    "toa_eff_temp",
    "toa_net_energy",
    "water_path",
    "wind_rot_div",
    "wind_speed",
    "zonal_mass_streamfunction",
)


def _precip_name_mapping(model=um):
    """Generate lists of variable names for `precip_sum()`."""
    return {
        "total": [model.ls_rain, model.ls_snow, model.cv_rain, model.cv_snow],
        "conv": [model.cv_rain, model.cv_snow],
        "stra": [model.ls_rain, model.ls_snow],
        "rain": [model.ls_rain, model.cv_rain],
        "snow": [model.ls_snow, model.cv_snow],
    }


@update_metadata(name="cos_lat", units="1")
def lat_cos(cube, model=um):
    """Convert the latitude coordinate to a cube and apply cosine to it."""
    lat_cube = coord_to_cube(cube, model.y, broadcast=False)
    lat_cos_cube = apply_ufunc(np.cos, apply_ufunc(np.deg2rad, lat_cube))
    return lat_cos_cube


@update_metadata(units="K")
def air_temperature(cubelist, const=None, model=um):
    """
    Get the real temperature from the given cube list.

    If not present, it is attempted to calculate it from other variables,
    Exner pressure or air pressure, and potential temperature.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes containing temperature.
    const: aeolus.const.const.ConstContainer, optional
        Must have `reference_surface_pressure` and `dry_air_gas_constant` as attributes.
        If not given, an attempt is made to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of air temperature.
    """
    try:
        return cubelist.extract_cube(model.temp)
    except ConMisErr:
        try:
            thta = cubelist.extract_cube(model.thta)
        except ConMisErr:
            raise MissingCubeError(f"Unable to get air temperature from {cubelist}")

        if len(cubelist.extract(model.exner)) == 1:
            exner = cubelist.extract_cube(model.exner)
        elif len(cubelist.extract(model.pres)) == 1:
            if const is None:
                const = thta.attributes["planet_conf"]
            pres = cubelist.extract_cube(model.pres)
            exner = (pres / const.reference_surface_pressure) ** (
                const.dry_air_gas_constant / const.dry_air_spec_heat_press
            ).data
        else:
            raise MissingCubeError(f"Unable to get air temperature from {cubelist}")
        temp = thta * exner
        temp.rename(model.temp)
        return temp


@update_metadata(units="K")
def air_potential_temperature(cubelist, const=None, model=um):
    """
    Get the potential temperature from the given cube list.

    If not present, it is attempted to calculate it from other variables,
    Exner pressure or air pressure, and real temperature.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes containing temperature.
    const: aeolus.const.const.ConstContainer, optional
        Must have `reference_surface_pressure` and `dry_air_gas_constant` as attributes.
        If not given, an attempt is made to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of air potential temperature.
    """
    try:
        return cubelist.extract_cube(model.thta)
    except ConMisErr:
        try:
            temp = cubelist.extract_cube(model.temp)
        except ConMisErr:
            raise MissingCubeError(f"Unable to get air potential temperature from {cubelist}")

        if len(cubelist.extract(model.exner)) == 1:
            exner = cubelist.extract_cube(model.exner)
        elif len(cubelist.extract(model.pres)) == 1:
            if const is None:
                const = temp.attributes["planet_conf"]
            pres = cubelist.extract_cube(model.pres)
            exner = (pres / const.reference_surface_pressure) ** (
                const.dry_air_gas_constant / const.dry_air_spec_heat_press
            ).data
        else:
            raise MissingCubeError(f"Unable to get air potential temperature from {cubelist}")
        thta = temp / exner
        thta.rename(model.thta)
        thta.convert_units("K")
        return thta


@const_from_attrs()
@update_metadata(units="Pa")
def air_pressure(cubelist, const=None, model=um):
    """
    Get pressure from the given cube list.

    If not present, it is attempted to calculate it from other variables,
    such as Exner pressure or air potential temperature and real temperature.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes containing related variables.
    const: aeolus.const.const.ConstContainer, optional
        Must have `reference_surface_pressure` and `dry_air_gas_constant` as attributes.
        If not given, an attempt is made to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of air pressure.
    """
    try:
        return cubelist.extract_cube(model.pres)
    except ConMisErr:
        try:
            exner = cubelist.extract_cube(model.exner)
            pres = (
                const.reference_surface_pressure
                * exner ** (const.dry_air_spec_heat_press / const.dry_air_gas_constant).data
            )
        except ConMisErr:
            try:
                temp = cubelist.extract_cube(model.temp)
                thta = cubelist.extract_cube(model.temp)
                exner = temp / thta
                pres = (
                    const.reference_surface_pressure
                    * exner ** (const.dry_air_spec_heat_press / const.dry_air_gas_constant).data
                )
            except ConMisErr:
                raise MissingCubeError(f"Unable to calculate air pressure from {cubelist}")
        pres.rename(model.pres)
        return pres


@update_metadata(units="kg m-3")
def air_density(cubelist, const=None, model=um):
    """
    Get air density from the given cube list.

    If not present, it is attempted to calculate it from pressure and temperature.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes containing temperature.
    const: aeolus.const.const.ConstContainer, optional
        Must have `dry_air_gas_constant` as an attribute.
        If not given, an attempt is made to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of air density.
    """
    try:
        return cubelist.extract_cube(model.dens)
    except ConMisErr:
        try:
            temp = cubelist.extract_cube(model.temp)
            pres = cubelist.extract_cube(model.pres)
            if const is None:
                const = pres.attributes["planet_conf"]
            rho = pres / (const.dry_air_gas_constant * temp)
            rho.rename(model.dens)
            return rho
        except ConMisErr:
            _msg = f"Unable to get variables from\n{cubelist}\nto calculate air density"
            raise MissingCubeError(_msg)


@const_from_attrs()
def calc_derived_cubes(cubelist, const=None, model=um):
    """Calculate additional variables."""
    try:
        cubelist.extract_cube(model.temp)
    except ConMisErr:
        cubelist.append(air_temperature(cubelist, const=const, model=model))
    try:
        cubelist.extract_cube(model.thta)
    except ConMisErr:
        cubelist.append(air_potential_temperature(cubelist, const=const, model=model))
    try:
        cubelist.extract_cube(model.pres)
    except ConMisErr:
        cubelist.append(air_pressure(cubelist, const=const, model=model))
    try:
        cubelist.extract_cube(model.dens)
    except ConMisErr:
        cubelist.append(air_density(cubelist, const=const, model=model))
    try:
        cubelist.extract_cube(model.ghgt)
    except ConMisErr:
        cubelist.append(geopotential_height(cubelist, const=const, model=model))


@update_metadata(units="m2 s-2")
def geopotential_height(cubelist, const=None, model=um):
    """
    Get geopotential height from the given cube list.

    If not present, the altitude coordinate is transformed into a cube.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes containing temperature.
    const: aeolus.const.const.ConstContainer, optional
        Must have `gravity` as an attribute.
        If not given, an attempt is made to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of geopotential height.
    """
    try:
        return cubelist.extract_cube(model.ghgt)
    except ConMisErr:
        try:
            cube_w_height = cubelist.extract(_dim_constr(model.z, strict=False))[0]
            if const is None:
                const = cube_w_height.attributes["planet_conf"]
            g_hgt = coord_to_cube(cube_w_height, model.z) * const.gravity
            g_hgt.attributes = {k: v for k, v in cube_w_height.attributes.items() if k != "STASH"}
            ensure_bounds(g_hgt, [model.z])
            g_hgt.rename(model.ghgt)
            return g_hgt
        except ConMisErr:
            _msg = f"No cubes in \n{cubelist}\nwith {model.z} as a coordinate."
            raise MissingCubeError(_msg)


def flux(cubelist, quantity, axis, weight_by_density=True, model=um):
    r"""
    Calculate horizontal or vertical flux of some quantity.

    .. math::
        F_{x} = u (\rho q),
        F_{y} = v (\rho q),
        F_{z} = w (\rho q)

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    quantity: str or iris.Constraint
        Quantity (present in the cube list).
    axis: str
        Axis of the flux component (x|y|z)
    weight_by_density: bool, optional
        Multiply by a cube of air density (must be present in the input cube list).
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of a flux component with the same dimensions as input cubes.
    """
    if axis.lower() == "x":
        u = cubelist.extract_cube(model.u)
    elif axis.lower() == "y":
        u = cubelist.extract_cube(model.v)
    elif axis.lower() == "z":
        u = cubelist.extract_cube(model.w)
    q = cubelist.extract_cube(quantity)
    fl = u * q
    if weight_by_density:
        fl *= cubelist.extract_cube(model.dens)
    fl.rename(f"flux_of_{quantity}_along_{axis.lower()}_axis")
    return fl


@update_metadata(units="W m-2")
def toa_cloud_radiative_effect(cubelist, kind, model=um):
    r"""
    Calculate domain-average TOA cloud radiative effect (CRE).

    .. math::
        CRE_{TOA} = F_{up,clear-sky} - F_{up,all-sky}

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes
    kind: str
        Shortwave ('sw'), longwave ('lw'), or 'total' CRE

    Returns
    -------
    iris.cube.Cube
        Cube of CRE.
    """
    name = f"toa_cloud_radiative_effect_{kind}"
    if kind == "sw":
        all_sky = model.toa_osr
        clr_sky = model.toa_osr_cs
    elif kind == "lw":
        all_sky = model.toa_olr
        clr_sky = model.toa_olr_cs
    elif kind == "total":
        sw = toa_cloud_radiative_effect(cubelist, "sw", model=model)
        lw = toa_cloud_radiative_effect(cubelist, "lw", model=model)
        cre = sw + lw
        cre.rename(name)
        return cre

    cube_clr = cubelist.extract_cube(clr_sky)
    cube_all = cubelist.extract_cube(all_sky)
    cre = cube_clr - cube_all
    cre.rename(name)
    return cre


@update_metadata(name="toa_net_downward_energy_flux", units="W m-2")
def toa_net_energy(cubelist, model=um):
    """
    Calculate domain-average TOA energy flux.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of total TOA downward energy flux.
    """
    varnames = [
        model.toa_isr,
        model.toa_osr,
        model.toa_olr,
    ]
    terms = cubelist.extract(varnames)
    if len(terms) != 3:
        raise MissingCubeError(
            f"{varnames} required for TOA energy balance are missing from cubelist:\n{cubelist}"
        )
    terms_ave = []
    for cube in terms:
        terms_ave.append(cube)
    toa_net = terms_ave[0] - terms_ave[1] - terms_ave[2]
    return toa_net


@update_metadata(name="surface_net_downward_energy_flux", units="W m-2")
def sfc_net_energy(cubelist, model=um):
    """
    Calculate domain-average surface energy flux.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes with net LW and SW radiation, sensible and latent surface fluxes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of total surface downward energy flux.
    """
    net_down_lw = cubelist.extract_cube(model.sfc_net_down_lw)
    net_down_sw = cubelist.extract_cube(model.sfc_net_down_sw)
    shf = cubelist.extract_cube(model.sfc_shf)
    lhf = cubelist.extract_cube(model.sfc_lhf)
    sfc_net = net_down_lw + net_down_sw - shf - lhf
    return sfc_net


@const_from_attrs()
@update_metadata(name="surface_net_downward_water_flux", units="mm h-1")
def sfc_water_balance(cubelist, const=None, model=um):
    """
    Calculate domain-average precipitation minus evaporation.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    const: aeolus.const.const.ConstContainer, optional
        Must have a scalar cube of `condensible_density` as an attribute.
        If not given, attempt is made to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of total surface downward water flux (P-E).
    """
    try:
        evap = cubelist.extract_cube(model.sfc_evap)
    except ConMisErr:
        try:
            lhf = cubelist.extract_cube(model.sfc_lhf)
            evap = lhf / const.condensible_heat_vaporization
            evap /= const.condensible_density
        except (KeyError, ConMisErr):
            raise MissingCubeError(f"Cannot retrieve evaporation from\n{cubelist}")
    try:
        precip = cubelist.extract_cube(model.ppn)
        precip /= const.condensible_density
    except ConMisErr:
        precip = precip_sum(cubelist, ptype="total", const=const, model=model)
    precip.convert_units("mm h-1")
    evap.convert_units("mm h-1")
    net = precip - evap
    return net


@const_from_attrs()
def precip_sum(cubelist, ptype="total", const=None, model=um):
    """
    Calculate a sum of different types of precipitation [:math:`mm~day^{-1}`].

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    ptype: str, optional
        Precipitation type (total|stra|conv|rain|snow).
    const: aeolus.const.const.ConstContainer, optional
        Must have a scalar cube of `condensible_density` as an attribute.
        If not given, attempt to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Sum of cubes of precipitation with units converted to mm per day.
    """
    try:
        varnames = _precip_name_mapping(model=model)[ptype]
    except KeyError:
        raise ArgumentError(f"Unknown ptype={ptype}")
    if len(cubelist.extract(varnames)) == 0:
        raise MissingCubeError(
            f"{varnames} required for ptype={ptype} are missing from cubelist:\n{cubelist}"
        )
    precip = 0.0
    for varname in varnames:
        try:
            cube = cubelist.extract_cube(varname)
            precip += cube
        except ConMisErr:
            pass
    precip /= const.condensible_density
    precip.convert_units("mm day^-1")
    precip.rename(f"{ptype}_precip_rate")
    return precip


@update_metadata(name="heat_redistribution_efficiency", units="1")
def heat_redist_eff(cubelist, region_a, region_b, model=um):
    r"""
    Heat redistribution efficiency (Leconte et al. 2013).

    .. math::
        \eta=\frac{OLR_{TOA,night}}{OLR_{TOA,day}}

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    region_a: aeolus.region.Region
        First region (usually nightside).
    region_b: aeolus.region.Region
        Second region (usually dayside).
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of eta parameter with collapsed spatial dimensions.
    """
    toa_olr = cubelist.extract_cube(model.toa_olr)
    toa_olr_a = spatial_mean(toa_olr.extract(region_a.constraint))
    toa_olr_b = spatial_mean(toa_olr.extract(region_b.constraint))
    eta = toa_olr_a / toa_olr_b
    return eta


@update_metadata(name="toa_effective_temperature", units="K")
def toa_eff_temp(cubelist, model=um):
    r"""
    Calculate effective temperature from TOA OLR.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of :math:`T_{eff}`.
    """
    toa_olr = cubelist.extract_cube(model.toa_olr)
    sbc = init_const().stefan_boltzmann
    t_eff = (toa_olr / sbc) ** 0.25
    return t_eff


@update_metadata(name="normalised_greenhouse_effect_parameter", units="1")
def ghe_norm(cubelist, model=um):
    r"""
    Normalised greenhouse effect parameter.

    .. math::
        GHE = 1 - \left(\frac{T_{eff}}{T_{sfc}}\right)^4

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of greenhouse effect parameter.
    """
    t_sfc = cubelist.extract_cube(model.t_sfc)
    t_eff = toa_eff_temp(cubelist, model=model)
    out = (t_eff / t_sfc) ** 4
    out = out.copy(data=1 - out.data)
    return out


@const_from_attrs()
@update_metadata(name="greenhouse_effect_parameter", units="K")
def greenhouse_effect(cubelist, kind="all_sky", const=None, model=um):
    r"""
    Calculate the greenhouse effect [K].

    .. math::
        GHE = T_{sfc} - \left(\frac{T_{eff}}{T_{sfc}}\right)^4

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    kind: str, optional
        Type of GHE:  "all_sky" or "clear_sky"
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of greenhouse effect parameter.
    """
    t_sfc = cubelist.extract_cube(model.t_sfc)
    if kind == "all_sky":
        toa_olr = cubelist.extract_cube(model.toa_olr)
    elif kind == "clear_sky":
        toa_olr = cubelist.extract_cube(model.toa_olr_cs)
    ghe = t_sfc - (toa_olr / const.stefan_boltzmann) ** 0.25
    return ghe


@const_from_attrs()
@update_metadata(name="bond_albedo", units="1")
def bond_albedo(cubelist, const=None, model=um):
    r"""
    Bold albedo.

    .. math::
        4 \frac{OSR_{TOA}}{S_{0}}

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    const: aeolus.const.const.ConstContainer, optional
        Must have a scalar cube of `solar_constant` as an attribute.
        If not given, attempt to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of bond albedo.
    """
    toa_osr = cubelist.extract_cube(model.toa_osr)
    sc = const.solar_constant
    b_alb = 4 * toa_osr / sc
    return b_alb


@update_metadata(units="kg m-2")
def water_path(cubelist, kind="water_vapour", model=um):
    r"""
    Water vapour or condensate path, i.e. a vertical integral of a water phase.

    .. math::
        WP = \int_{z_{sfc}}^{z_{top}} \rho q dz

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes containing appropriate mixing ratio and air density.
    kind: str, optional
        Short name of the water phase to be integrated.
        Options are water_vapour (default) | liquid_water | ice_water | cloud_water
        `cloud_water` is the sum of liquid and ice phases.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.
        `model.z` is used as a vertical coordinate for integration.

    Returns
    -------
    iris.cube.Cube
        Cube of water path with collapsed vertical dimension.
    """
    if kind == "water_vapour":
        q = cubelist.extract_cube(model.sh)
    elif kind == "liquid_water":
        q = cubelist.extract_cube(model.cld_liq_mf)
    elif kind == "ice_water":
        q = cubelist.extract_cube(model.cld_ice_mf)
    elif kind == "cloud_water":
        q = cubelist.extract_cube(model.cld_liq_mf).copy()
        q += cubelist.extract_cube(model.cld_ice_mf)
    rho = cubelist.extract_cube(model.dens)
    wp = integrate(q * rho, model.z)
    wp.rename(f"{kind}_path")
    return wp


def dry_lapse_rate(cubelist, model=um):
    r"""
    Dry lapse rate, or the change of air temperature with altitude.

    .. math::
        \gamma = \partial T_{air} / \partial z

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes containing temperature.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of dT/dz.
    """
    temp = cubelist.extract_cube(model.temp)
    return d_dz(temp, model=model)


def horiz_wind_cmpnts(cubelist, model=um):
    """
    Extract u and v wind components from a cube list and interpolate v on u's grid if necessary.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        List of cubes with horizontal wind components
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    u, v: iris.cube.Cube
        Cubes of wind components.
    """
    # use relaxed extracting
    u = cubelist.extract(model.u)[0]
    v = cubelist.extract(model.v)[0]
    # interpolate v on u's grid if coordinates are different
    v = regrid_3d(v, u, model=model)
    return u, v


@const_from_attrs()
@update_metadata(name="local_superrotation_index", units="1")
def superrotation_index(cubelist, const=None, model=um):
    r"""
    Local superrotation index.

    .. math::
        s = \frac{m}{\Omega a^2} - 1,

        m = a cos\phi(\Omega a cos\phi + u)

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        List of cubes containing a cube of zonal velocity (u).
    const: aeolus.const.const.ConstContainer, optional
        Constainer with the relevant planetary constants.
        If not given, attempt is made to retrieve it from cube attributes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    s_idx: iris.cube.Cube
        Cubes of superrotation index.

    References
    ----------
    Read (1986), https://doi.org/10.1002/qj.49711247114
    """
    # Zonal velocity
    u = cubelist.extract_cube(model.u)
    # Radius of the planet
    r = const.radius
    r.convert_units("m")
    # Rotation rate
    omega = const.planet_rotation_rate

    lat_coord = u.coord(model.y).copy()
    lat_dim = u.coord_dims(model.y)[0]
    lat_coord.convert_units("radians")
    lat_cos_coord = _coord_cos(lat_coord)

    # Calculate intermediate terms
    omega_r_coslat = omega.data * r.data * lat_cos_coord
    omega_r_coslat.units = Unit("m s-1")
    r_coslat = r.data * lat_cos_coord
    r_coslat.units = Unit("m")
    inner_sum = add(u, omega_r_coslat, dim=lat_dim)

    # Calculate axial component of specific absolute angular momentum
    ang_mom = multiply(inner_sum, r_coslat, dim=lat_dim)

    # Final index
    s_idx = ang_mom / omega / (r**2)
    s_idx.convert_units("1")
    s_idx = s_idx.copy(data=s_idx.data - 1)
    return s_idx


@update_metadata(name="wind_speed", units="m s-1")
def wind_speed(*components):
    r"""
    Calculate the wind speed (magnitude of the wind vector).

    .. math::
        \sqrt{u^2 + v^2 + w^2}

    Parameters
    ----------
    args: iris.cube.Cube
        Cubes of u, v, w wind components.

    Returns
    -------
    iris.cube.Cube
    """
    out = sum(cube**2 for cube in components) ** 0.5
    return out


@update_metadata(name="atmosphere_hybrid_sigma_pressure_coordinate", units="1")
def sigma_p(cubelist, model=um):
    r"""
    Calculate sigma (normalised pressure coordinate) from a cube of pressure.

    .. math::
        \sigma = p / p_{sfc}

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.
    """
    pres_atm = cubelist.extract_cube(model.pres)
    pres_sfc = cubelist.extract_cube(model.p_sfc)
    out = pres_atm / pres_sfc
    return out


@const_from_attrs()
@update_metadata(name="zonal_mass_streamfunction", units="kg s^-1")
def zonal_mass_streamfunction(cubelist, const=None, model=um):
    r"""
    Calculate mean zonal mass streamfunction.

    .. math::
        \Psi_Z = 2\pi a \int_{z_{sfc}}^{z_{top}}\overline{\rho}^* \overline{u}^* dz

    References
    ----------
    Haqq-Misra & Kopparapu (2015), eq. 5;
    Hartmann (1994), Global Physical Climatology, eq. 6.21

    Examples
    --------
    >>> lat_band_constr = iris.Constraint(
        latitude=lambda x: -30 <= x.point <= 30, longitude=lambda x: True
    )
    >>> mzsf = zonal_mass_streamfunction(cubes.extract(lat_band_constr))
    """
    streamf_const = 2 * np.pi * const.radius

    u = cubelist.extract_cube(model.u)
    u_tm = time_mean(u, model=model)
    u = -1 * (u_tm - zonal_mean(u_tm, model=model))

    if u.coord(axis="z").units.is_convertible("m"):
        rho = cubelist.extract_cube(model.dens).copy()
        rho = time_mean(rho, model=model)
        rho.coord(model.z).bounds = None
        u.coord(model.z).bounds = None
        integrand = u * rho
        # integrand = reverse(integrand, model.z)
        # print(integrand.coord(um.z))
        res = cumsum(integrand, "z", axis_weights=True, model=model)
        # res = cumsum(integrand, "z", axis_weights=True, model=model)
        # res = reverse(res, model.z)

    elif u.coord(axis="z").units.is_convertible("Pa"):
        integrand = u
        res = cumsum(integrand, "z", axis_weights=True, model=model)
        res /= const.gravity

    res *= streamf_const
    return res


@const_from_attrs()
@update_metadata(name="meridional_mass_streamfunction", units="kg s^-1")
def meridional_mass_streamfunction(cubelist, const=None, model=um):
    r"""
    Calculate the mean meridional mass streamfunction.

    * In height coordinates

    .. math::
        \Psi_M = - 2\pi cos\phi a \int_{z_{sfc}}^{z_{top}}\overline{\rho v} dz

    * In pressure coordinates

    .. math::
        \Psi_M = 2\pi cos\phi a \int_{0}^{p_{sfc}}\overline{v} dp / g

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input cubelist.
    const: aeolus.const.const.ConstContainer, optional
        If not given, constants are attempted to be retrieved from
        attributes of a cube in the cube list.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube with collapsed spatial dimensions.

    References
    ----------
    Haqq-Misra & Kopparapu (2015), eq. 4;
    Vallis (2017)

    Examples
    --------
    >>> from aeolus.calc import meridional_mass_streamfunction, time_mean
    >>> from aeolus.const import init_const
    >>> from aeolus.model import um
    >>> earth_constants = init_const("earth")
    >>> cubes = iris.cube.CubeList([time_mean(cube) for cube in input_cubes])
    >>> mmsf = meridional_mass_streamfunction(cubes, const=earth_constants, model=um)
    """
    v = cubelist.extract_cube(model.v)
    v = zonal_mean(v, model=model)
    if v.coord(model.z).units.is_convertible("m"):
        rho = zonal_mean(cubelist.extract_cube(model.dens), model=model)
        rho.coord(model.z).bounds = None
        v.coord(model.z).bounds = None
        # Reverse the coordinate to start from the model top (where p=0 or z=z_top)
        # TODO: check if the coordinate is ascending or descending
        integrand = reverse(v * rho, model.z)
        res = -1 * cumsum(integrand, "z", axis_weights=True, model=model)
        # Reverse the result back
        res = reverse(res, model.z)
    elif v.coord(model.z).units.is_convertible("Pa"):
        res = cumsum(v, "z", axis_weights=True, model=model)
        res /= const.gravity
    # Calculate the constant: 2 pi cos\phi a
    streamf_const = 2 * np.pi * const.radius * lat_cos(res, model=model)
    res *= streamf_const
    return res


@const_from_attrs()
def wind_rot_div(u, v, truncation=None, const=None, model=um):
    """
    Split the horizontal wind field into divergent and rotational parts (Helmholtz decomposition).

    The Helmholtz decomposition method uses the `windspharm` library:
    https://ajdawson.github.io/windspharm/latest/

    Parameters
    ----------
    u: iris.cube.Cube
        Eastward wind.
    v: iris.cube.Cube
        Northward wind.
    truncation: int
        Truncation for the spherical harmonic computation (See windspharm docs for details).
    const: aeolus.const.const.ConstContainer, optional
        If not given, constants are attempted to be retrieved from
        attributes of a cube in the cube list.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    out: dict
        Dictionary of cubes of:
          - input wind components (for convenience),
          - divergent components,
          - rotational components,
          - zonal mean rotational components,
          - zonal eddy rotational components.

    References
    ----------
    Hammond and Lewis (2021), https://doi.org/10.1073/pnas.2022705118.
    """
    from windspharm.iris import VectorWind

    vec = VectorWind(u, v, rsphere=const.radius.data)
    div_cmpnt_u, div_cmpnt_v, rot_cmpnt_u, rot_cmpnt_v = vec.helmholtz(truncation=truncation)
    out = {}
    out["u_total"] = u
    out["v_total"] = v
    out["u_div"] = reverse(div_cmpnt_u, model.y)
    out["v_div"] = reverse(div_cmpnt_v, model.y)
    out["u_rot"] = reverse(rot_cmpnt_u, model.y)
    out["v_rot"] = reverse(rot_cmpnt_v, model.y)

    for cmpnt in ["u", "v"]:
        rot_cmpnt = out[f"{cmpnt}_rot"]
        out[f"{cmpnt}_rot_zm"] = preserve_shape(zonal_mean)(rot_cmpnt)
        out[f"{cmpnt}_rot_zm"].rename(f"zonal_mean_of_{rot_cmpnt.name()}")
        out[f"{cmpnt}_rot_eddy"] = rot_cmpnt - out[f"{cmpnt}_rot_zm"]
        out[f"{cmpnt}_rot_eddy"].rename(f"zonal_deviation_of_{rot_cmpnt.name()}")
    return out
