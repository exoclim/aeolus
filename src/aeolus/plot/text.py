# -*- coding: utf-8 -*-
"""Text-related formatting functions."""
import itertools
import math
import string

from ..exceptions import ArgumentError
from ..calc import spatial_mean


__all__ = (
    "cube_minmeanmax_str",
    "fmt_lonlat",
    "subplot_label_generator",
    "tex2cf_units",
    "unit_format",
)


def cube_minmeanmax_str(
    cube,
    sep=" | ",
    eq_sign="=",
    fmt="auto",
    weight=True,
    labels=["min", "mean", "max"],
    **kw_unit_format,
):
    """Return min, mean and max of an iris cube as a string."""
    # Compute the stats
    _min = float(cube.data.min())
    if weight:
        _mean = float(spatial_mean(cube).data)
    else:
        _mean = float(cube.data.mean())
    _max = float(cube.data.max())
    # Assemble a string
    txts = []
    for label, num in zip(labels, [_min, _mean, _max]):
        if fmt == "auto":
            if (math.log10(abs(_mean)) < 0) or (math.log10(abs(_mean)) > 5):
                _str = f"{label}{eq_sign}{num:.0e}"
            else:
                _str = f"{label}{eq_sign}{round(num):.0f}"
        elif fmt == "pretty":
            _str = f"{label}{eq_sign}{unit_format(num, **kw_unit_format)}"
        else:
            _str = f"{label}{eq_sign}{num:{fmt}}"
        txts.append(_str)
    return sep.join(txts)


def fmt_lonlat(value, lon_or_lat, degree=False):
    r"""
    Convert longitude or latitude value to string with a hemisphere identifier.

    Parameters
    ----------
    value: int
        Value of longitude or latitude. Note that this function is only for integer values.
    lon_or_lat: str
        Longitude or latitude
    degree: bool, optional
        If true, a TeX degree symbol is included

    Returns
    -------
    str

    Examples
    --------
    >>> fmt_lonlat(-25, "lon")
    '25W'
    >>> fmt_lonlat(89, "lat", degree=True)
    '89$^\\degree$N'
    >>> fmt_lonlat(0, "lon")
    '0'
    """
    from LatLon23 import Latitude, Longitude  # noqa

    if lon_or_lat.lower().startswith("lat"):
        res = Latitude(value)
    elif lon_or_lat.lower().startswith("lon"):
        res = Longitude(value)
    else:
        raise ArgumentError("2nd arg or the function should start with `lon` or `lat`")
    out = res.to_string("%d%H")
    if degree:
        out = out[:-1] + r"$^\degree$" + out[-1]
    if value == 0:
        out = out[:-1]
    return out


def subplot_label_generator():
    """Return generator of alphabetic labelling of subplots."""
    for i in itertools.count(1):
        for p in itertools.product(string.ascii_lowercase, repeat=i):
            yield "".join(p)


def tex2cf_units(unit_str):
    """Convert a TeX string to a string that can be used in cf_units."""
    return unit_str.replace("$", "").replace("{", "").replace("}", "").replace("^", "**")


def unit_format(value, unit="1", decimal_digits=1, precision=None, exponent=None):
    r"""
    Return a string representation of a given number with units.

    Format the scientific notation of the given number for use with LaTeX or Mathtext,
    with specified number of significant decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified explicitly.

    Parameters
    ----------
    value: int or float
        Number to be formatted.
    unit: str or int, optional
        Units of the number, if any.
    decimal_digits: int, optional
        Significant decimal digits.
    precision: int, optional
        Number of decimal digits to show.
    exponent: int, optional
        Exponent of the resulting number.

    Returns
    -------
    str

    Examples
    --------
    >>> unit_format(-1.234e-5)
    '$-1.2\\times10^{-5}$'
    >>> unit_format(6)
    '$6.0$'
    >>> unit_format(1.234e5, unit="m s^{-1}")
    '$1.2\\times10^{5}$ $m$ $s^{-1}$'
    >>> unit_format(1.234, decimal_digits=1, precision=3)
    '$1.200$'
    >>> unit_format(1.234, decimal_digits=3, precision=1)
    '$1.2$'
    """
    if unit in ["", "1", 1]:
        unit_str = ""
    else:
        unit_str = rf'${str(unit).replace(" ", "$ $")}$'
        if "%" in unit_str:
            unit_str = unit_str.replace("%", r"\%")
    if value == 1:
        string = unit_str
    else:
        if precision is None:
            precision = decimal_digits
        if exponent is None:
            exponent = int(math.floor(math.log10(abs(value))))

        coeff = round(value / 10**exponent, decimal_digits)

        if exponent in [0, 1]:
            string = rf"${round(value, decimal_digits):.{precision}f}$"
        else:
            string = rf"${coeff:.{precision}f}\times10^{{{exponent:d}}}$"
        if not unit == "1":
            string += rf" {unit_str}"
    return string
