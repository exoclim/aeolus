"""Functions to calculate commonly used variables in atmospheric science."""
from warnings import warn

import iris

import numpy as np

from .coord_utils import UM_LATLON, UM_TIME, ensure_bounds
from .exceptions import AeolusWarning, ArgumentError
from .grid import area_weights_cube
from .subset import extract_last_year


__all__ = (
    "calc_spatial",
    "calc_spatial_quartiles",
    "calc_meridional_mean",
    "calc_toa_cloud_radiative_effect",
    "calc_toa_net",
    "calc_total_precip",
    "calc_sfc_water_balance",
    "integrate",
)


def calc_spatial(cube, aggr, coords=UM_LATLON):
    """
    Calculate spatial statistic with geographic grid weights.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with longitude and latitude coordinates.
    aggr: str
        Statistical aggregator (see iris.analysis for available aggregators).
    coords: list, optional
        List of names of spatial coordinates.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.

    Examples
    --------
    >>> calc_spatial(my_data_cube, "mean")

    """
    ensure_bounds(cube)
    flag = all(cube.coord(c).has_bounds() for c in coords)
    aggregator = getattr(iris.analysis, aggr.upper())
    if flag and isinstance(aggregator, iris.analysis.WeightedAggregator):
        kw = {"weights": area_weights_cube(cube, normalize=True).data}
    else:
        kw = {}
    return cube.collapsed(coords, aggregator, **kw)


def calc_spatial_quartiles(cube):
    """Calculate quartiles over horizontal coordinates."""
    warn("No weights are applied!", AeolusWarning)
    q25 = cube.collapsed(UM_LATLON, iris.analysis.PERCENTILE, percent=25)
    q75 = cube.collapsed(UM_LATLON, iris.analysis.PERCENTILE, percent=75)
    return q25, q75


def calc_meridional_mean(cube, lat_name=UM_LATLON[0]):
    """
    Calculate cube's meridional average.

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube with a latitude coordinate.
    lat_name: str, optional
        Name of the latitude coordinate.

    Returns
    -------
    iris.cube.Cube
        Collapsed cube.
    """
    coslat = np.cos(np.deg2rad(cube.coord(lat_name).points))
    coslat2d = iris.util.broadcast_to_shape(coslat, cube.shape, (0,))
    cube_mean = (cube * coslat2d).collapsed(lat_name, iris.analysis.SUM) / np.sum(coslat)
    return cube_mean


def last_year_mean(cube):
    """Get the time mean of over the last year."""
    last_year_constraint = extract_last_year(cube)
    return cube.extract(last_year_constraint).collapsed(UM_TIME, iris.analysis.MEAN)


def calc_toa_cloud_radiative_effect(cubelist, kind):
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
        sw = calc_toa_cloud_radiative_effect(cubelist, "sw")
        lw = calc_toa_cloud_radiative_effect(cubelist, "lw")
        cre = sw + lw
        cre.rename(name)
        return cre

    cube_clr = calc_spatial(
        cubelist.extract_strict(iris.AttributeConstraint(STASH=clr_sky)), "mean"
    )
    cube_all = calc_spatial(
        cubelist.extract_strict(iris.AttributeConstraint(STASH=all_sky)), "mean"
    )
    cre = cube_clr - cube_all
    cre.rename(name)
    return cre


def calc_toa_net(cubelist):
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
        terms_ave.append(calc_spatial(cube, "mean"))
    toa_net = terms_ave[0] - terms_ave[1] - terms_ave[2]
    toa_net.rename("toa_net_radiative_flux")
    return toa_net


def calc_sfc_water_balance(cubelist):
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
        terms_ave.append(calc_spatial(cube, "mean"))
    net = terms_ave[0] - terms_ave[1]
    net.rename("surface_net_water_flux")
    return net


def calc_total_precip(cubelist, ptype=None):
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


def integrate(cube, coord):
    """
    Integrate the cube along a 1D coordinate using the trapezoidal rule.

    Note: `coord` must be one of the dimensional coordinates of the cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube containing the given coordinate.
    coord: str or iris.coords.Coord
        Coordinate for integration.

    Returns
    -------
    iris.cube.Cube
        integrated cube.
    """
    # TODO: allow non-dim coordinates
    c = cube.coord(coord)
    others = [dc.name() for dc in cube.dim_coords if dc.name() != c.name()]
    dim = cube.coord_dims(c)[0]
    data = np.trapz(cube.data, c.points, axis=dim)
    res = next(cube.slices(others)).copy(data=data)
    res.units = cube.units * c.units
    res.remove_coord(c)
    res.rename(f"integral_of_{cube.name()}_wrt_{c.name()}")
    # ensure_bounds(cube, [c])
    # delta = iris.coords.AuxCoord(c.bounds[:, 1] - c.bounds[:, 0], units=c.units)
    # res = iris.analysis.maths.multiply(cube, delta, dim=dim).collapsed(c, iris.analysis.SUM)
    return res
