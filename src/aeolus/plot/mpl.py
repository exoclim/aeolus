# -*- coding: utf-8 -*-
"""Matplotlib-related utilities."""
from typing import Sequence, Union, Optional

from iris.cube import Cube
from matplotlib.axes._axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

from ..model import lfric, um
from ..model.base import Model
from ..coord import get_cube_rel_days

__all__ = (
    "add_custom_legend",
    "capitalise",
    "hcross",
    "make_list_2d",
    "map_scatter",
    "timeseries_1d",
    "timeseries_2d",
)


def add_custom_legend(ax_or_fig: Union[Axes, Figure], styles_and_labels: dict, **leg_kw) -> None:
    """
    Add a custom legend to a matplotlib axis or figure.

    Parameters
    ----------
    ax_or_fig: matplotlib.axes._subplots.AxesSubplot or matplotlib.figure.Figure
        Matplotlib object where to put the legend.
    styles_and_labels: dict
        Dictionary with labels as keys and a dictionary of plot
        keywords as values.
    leg_kw: dict, optional
        Keyword arguments passed to `legend()` function.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>> ax = plt.axes()
    >>> my_dict = dict(foo=dict(color='C0', marker="X"),
                       bar=dict(color='C1', marker="o"))
    >>> add_custom_legend(ax, my_dict, loc=2, title="blah")

    """
    lines = [Line2D([0], [0], **style) for style in styles_and_labels.values()]
    leg = ax_or_fig.legend(lines, styles_and_labels.keys(), **leg_kw)
    try:
        if ax_or_fig.legend_ is not None:
            ax_or_fig.add_artist(leg)
    except AttributeError:
        pass


def capitalise(s: str, sep_old: Optional[str] = "_", sep_new: Optional[str] = " ") -> str:
    """Split the string and capitalise each word."""
    return sep_new.join([i.capitalize() for i in s.split(sep_old)])


def hcross(
    cube: Cube, ax: Optional[Axes] = None, model: Optional[Model] = um, **kw_plt
) -> Optional[Axes]:
    """Plot a horizontal cross-section aka lat-lon map of a 2D cube."""
    newax = False
    if ax is None:
        ax = plt.axes()
        newax = True
    fig = ax.figure
    lons = cube.coord(um.x).points
    lats = cube.coord(um.y).points
    mappable = ax.pcolormesh(lons, lats, cube.data, **kw_plt)
    fig.colorbar(mappable, ax=ax)
    if newax:
        return ax


def make_list_2d(
    list_x: Sequence[str], list_y: Sequence[str], transpose: Optional[bool] = False
) -> list:
    """Create a nested list out of 2 given lists."""
    if transpose:
        return [[f"{key_y}-{key_x}" for key_x in list_x] for key_y in list_y]
    else:
        return [[f"{key_x}-{key_y}" for key_x in list_x] for key_y in list_y]


def map_scatter(
    cube: Cube, ax: Optional[Axes] = None, model: Optional[Model] = lfric, **kw_plt
) -> Optional[Axes]:
    """Plot a lat-lon scatter plot of a 2D cube."""
    newax = False
    if ax is None:
        ax = plt.axes()
        newax = True
    fig = ax.figure
    # This doesn't work because lons and lats in the mesh are of size N+2
    # while the data array is of size N
    # lons, lats = cube.mesh.node_coords
    # lons, lats = lons.points, lats.points
    lons = cube.coord(model.x).points
    lats = cube.coord(model.y).points
    mappable = ax.scatter(lons, lats, c=cube.data, **kw_plt)
    fig.colorbar(mappable, ax=ax)
    if newax:
        return ax


def timeseries_1d(
    cube: Cube, ax: Optional[Axes] = None, model: Optional[Model] = um, **kw_plt
) -> Optional[Axes]:
    """Plot time series of a 1D cube."""
    newax = False
    if ax is None:
        ax = plt.axes()
        newax = True
    days = get_cube_rel_days(cube, model=um)
    ax.plot(days, cube.data, **kw_plt)
    if newax:
        return ax


def timeseries_2d(
    cube: Cube, ax: Optional[Axes] = None, model: Optional[Model] = um, **kw_plt
) -> Optional[Axes]:
    """Plot time series of a 2D cube."""
    newax = False
    if ax is None:
        ax = plt.axes()
        newax = True
    fig = ax.figure
    days = get_cube_rel_days(cube, model=um)
    z = cube.coord(um.z).points
    mappable = ax.pcolormesh(days, z, cube.data.T, **kw_plt)
    fig.colorbar(mappable, ax=ax)
    if newax:
        return ax
