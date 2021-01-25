"""Some commonly used diagnostics in atmospheric science."""
from cf_units import Unit

from iris.analysis.calculus import _coord_cos
from iris.analysis.maths import add, multiply
from iris.exceptions import ConstraintMismatchError as ConMisErr

import numpy as np

from .calculus import d_dz, integrate
from .meta import const_from_attrs, update_metadata
from .stats import spatial_mean
from ..const import init_const
from ..const.const import ScalarCube
from ..coord import coord_to_cube, ensure_bounds, regrid_3d
from ..exceptions import ArgumentError, MissingCubeError
from ..model import um
from ..subset import _dim_constr


__all__ = (
    "air_density",
    "air_potential_temperature",
    "air_temperature",
    "bond_albedo",
    "dry_lapse_rate",
    "flux",
    "geopotential_height",
    "ghe_norm",
    "heat_redist_eff",
    "horiz_wind_cmpnts",
    "precip_sum",
    "sfc_net_energy",
    "sfc_water_balance",
    "toa_cloud_radiative_effect",
    "toa_eff_temp",
    "toa_net_energy",
    "water_path",
    "wind_speed",
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
            exner = (pres / const.reference_surface_pressure.asc) ** (
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
            exner = (pres / const.reference_surface_pressure.asc) ** (
                const.dry_air_gas_constant / const.dry_air_spec_heat_press
            ).data
        else:
            raise MissingCubeError(f"Unable to get air potential temperature from {cubelist}")
        thta = temp / exner
        thta.rename(model.thta)
        thta.convert_units("K")
        return thta


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
            rho = pres / (const.dry_air_gas_constant.asc * temp)
            rho.rename(model.dens)
            return rho
        except ConMisErr:
            _msg = f"Unable to get variables from\n{cubelist}\nto calculate air density"
            raise MissingCubeError(_msg)


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
            g_hgt = coord_to_cube(cube_w_height, model.z) * const.gravity.asc
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
        sw = toa_cloud_radiative_effect(cubelist, "sw")
        lw = toa_cloud_radiative_effect(cubelist, "lw")
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


@const_from_attrs
@update_metadata(name="surface_net_downward_water_flux", units="mm h-1")
def sfc_water_balance(cubelist, const=None, model=um):
    """
    Calculate domain-average precipitation minus evaporation.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    const: aeolus.const.const.ConstContainer, optional
        Must have a `ScalarCube` of `condensible_density` as an attribute.
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
            evap = lhf / const.condensible_heat_vaporization.asc
            evap /= const.condensible_density.asc
        except (KeyError, ConMisErr):
            raise MissingCubeError(f"Cannot retrieve evaporation from\n{cubelist}")
    try:
        precip = cubelist.extract_cube(model.ppn)
        precip /= const.condensible_density.asc
    except ConMisErr:
        precip = precip_sum(cubelist, ptype="total", const=const, model=model)
    precip.convert_units("mm h-1")
    evap.convert_units("mm h-1")
    net = precip - evap
    return net


@const_from_attrs
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
        Must have a `ScalarCube` of `condensible_density` as an attribute.
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
    precip /= const.condensible_density.asc
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
    try:
        t_eff = (toa_olr / sbc) ** 0.25
    except ValueError:
        t_eff = (toa_olr / sbc.asc) ** 0.25
    return t_eff


@update_metadata(name="normalised_greenhouse_effect_parameter", units="1")
def ghe_norm(cubelist, model=um):
    r"""
    Normalised greenhouse effect parameter.

    .. math::
        GHE = 1 - \left(\frac{T_{eff}}{T_{sfc}}\right)^{1/4}

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        Cube of greenhouse effect parameter with collapsed spatial dimensions.
    """
    t_sfc = spatial_mean(cubelist.extract_cube(um.t_sfc))
    t_eff = toa_eff_temp(cubelist, model=model)
    one = t_eff.copy(data=np.ones(t_eff.shape))
    one.units = "1"
    gh_norm = one - (t_eff / t_sfc) ** 4
    return gh_norm


@const_from_attrs
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
        Must have a `ScalarCube` of `condensible_density` as an attribute.
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
    try:
        b_alb = 4 * toa_osr / sc
    except ValueError:
        b_alb = 4 * toa_osr / sc.asc
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


@const_from_attrs
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
    omega = ScalarCube.from_cube((const.day / (2 * np.pi)) ** (-1))
    omega.convert_units("s-1")

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
    s_idx = ang_mom / omega.asc / r.asc / r.asc
    s_idx.convert_units("1")
    s_idx = s_idx.copy(data=s_idx.data - 1)
    return s_idx


@update_metadata(name="wind_speed", units="m s-1")
def wind_speed(*components):
    r"""
    Calculate the wind speed (magnitude of the wind vector).

    Parameters
    ----------
    args: iris.cube.Cube
        Cubes of u, v, w wind components.

    .. math::
        \sqrt{u^2 + v^2 + w^2}
    """
    out = sum(cube ** 2 for cube in components) ** 0.5
    return out
