"""Some commonly used diagnostics in atmospheric science."""
import iris

from .stats import spatial
from ..exceptions import ArgumentError


__all__ = ("toa_cloud_radiative_effect", "toa_net", "total_precip", "sfc_water_balance")


def toa_cloud_radiative_effect(cubelist, kind):
    """
    Calculate domain-average TOA cloud radiative effect (CRE).

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes
    kind: str
        Shortwave ('sw'), longwave ('lw'), or 'total' CRE

    Returns
    -------
    iris.cube.Cube
        Cube of CRE with reduced dimensions.
    """
    name = f"toa_cloud_radiative_effect_{kind}"
    if kind == "sw":
        all_sky = "m01s01i208"
        clr_sky = "m01s01i209"
    elif kind == "lw":
        all_sky = "m01s02i205"
        clr_sky = "m01s02i206"
    elif kind == "total":
        sw = toa_cloud_radiative_effect(cubelist, "sw")
        lw = toa_cloud_radiative_effect(cubelist, "lw")
        cre = sw + lw
        cre.rename(name)
        return cre

    cube_clr = spatial(cubelist.extract_strict(iris.AttributeConstraint(STASH=clr_sky)), "mean")
    cube_all = spatial(cubelist.extract_strict(iris.AttributeConstraint(STASH=all_sky)), "mean")
    cre = cube_clr - cube_all
    cre.rename(name)
    return cre


def toa_net(cubelist):
    """
    Calculate domain-average TOA radiative flux.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes

    Returns
    -------
    iris.cube.Cube
        Cube with reduced dimensions.
    """
    terms = cubelist.extract(
        ["toa_incoming_shortwave_flux", "toa_outgoing_shortwave_flux", "toa_outgoing_longwave_flux"]
    )
    assert len(terms) == 3, f"Error when extracting TOA radiation terms from\n{cubelist}"
    terms_ave = []
    for cube in terms:
        terms_ave.append(spatial(cube, "mean"))
    toa_net = terms_ave[0] - terms_ave[1] - terms_ave[2]
    toa_net.rename("toa_net_radiative_flux")
    return toa_net


def sfc_water_balance(cubelist):
    """
    Calculate domain-average precipitation minus evaporation.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes

    Returns
    -------
    iris.cube.Cube
        Cube with reduced dimensions.
    """
    terms = cubelist.extract(["precipitation_flux", "surface_upward_water_flux"])
    terms_ave = []
    for cube in terms:
        terms_ave.append(spatial(cube, "mean"))
    net = terms_ave[0] - terms_ave[1]
    net.rename("surface_net_water_flux")
    return net


def total_precip(cubelist, ptype=None):
    """
    Calculate total precipitation flux [mm day^-1].

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Input list of cubes.
    ptype: str, optional
        Precipitation type (stra|conv).

    Returns
    -------
    iris.cube.Cube
        Sum of cubes of precipitation with units converted to mm per day.
    """
    conv = ["convective_rainfall_flux", "convective_snowfall_flux"]
    stra = ["stratiform_rainfall_flux", "stratiform_snowfall_flux"]
    if ptype is None:
        varnames = stra + conv
    elif ptype == "stra":
        varnames = stra
    elif ptype == "conv":
        varnames = conv
    else:
        raise ArgumentError(f"Unknown ptype={ptype}")
    tot_precip = cubelist.extract_strict(varnames[0]).copy()
    for varname in varnames[1:]:
        try:
            tot_precip += cubelist.extract_strict(varname)
        except iris.exceptions.ConstraintMismatchError:
            pass
    tot_precip /= iris.coords.AuxCoord(1000, units="kg m-3")
    tot_precip.convert_units("mm day-1")
    tot_precip.rename("total_precipitation_flux")
    return tot_precip
