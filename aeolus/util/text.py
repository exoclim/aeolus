"""Text-related formatting functions."""
import LatLon23 as LL

from ..exceptions import ArgumentError


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
    if lon_or_lat.lower().startswith("lat"):
        res = LL.Latitude(value)
    elif lon_or_lat.lower().startswith("lon"):
        res = LL.Longitude(value)
    else:
        raise ArgumentError("2nd arg or the function should start with `lon` or `lat`")
    out = res.to_string("%d%H")
    if degree:
        out = out[:-1] + r"$^\degree$" + out[-1]
    if value == 0:
        out = out[:-1]
    return out
