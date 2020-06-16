"""Plotting functions used with cartopy."""
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes

from matplotlib.transforms import offset_copy

from mpl_toolkits.axes_grid1 import AxesGrid

from .text import fmt_lonlat


__all__ = ("GeoAxesGrid", "label_global_map_gridlines")


class GeoAxesGrid(AxesGrid):
    """
    Grid of cartopy axes.

    A subclass of :class:`mpl_toolkits.axes_grid1.AxesGrid` representing
    a grid of maps with the same projection :class:`~cartopy.crs.Projection`.

    - `axes_class` is defined automatically
    - The :class:`AxesGrid` built-in labelling is always switched off,
      and instead a standard procedure of creating
      grid lines and labels should be used.
    """

    def __init__(self, fig, rect, nrows_ncols, projection, **axesgrid_kw):
        """
        Initialise GeoAxesGrid.

        Build a :class:`GeoAxesGrid` instance with a grid nrows*ncols
        :class:`GeoAxes` with a projection :class:`~cartopy.crs.Projection`
        in :class:`~matplotlib.figure.Figure` *fig* with
        *rect=[left, bottom, width, height]* (in
        :class:`~matplotlib.figure.Figure` coordinates) or
        the subplot position code (e.g., "121").
        """
        axesgrid_kw["axes_class"] = (GeoAxes, {"map_projection": projection})
        axesgrid_kw["label_mode"] = ""  # note the empty label_mode
        super(GeoAxesGrid, self).__init__(fig, rect, nrows_ncols, **axesgrid_kw)


def label_global_map_gridlines(
    fig, ax, xticks=[], yticks=[], xoff=-10, yoff=-10, degree=False, **text_kw
):
    """
    Label gridlines of a global cartopy map.

    Parameters
    ----------
    fig: matplotlib.figure.Figure
        Figure object.
    ax: cartopy.mpl.geoaxes.GeoAxesSubplot
        Cartopy axes.
    xticks: array-like, optional
        Sequence of longitude ticks.
    yticks: array-like, optional
        Sequence of latitude ticks.
    xoff: float, optional
        Longitude label offset from the axis (units are points).
        If negative (by default), the labels are drawn at the east boundary,
        otherwise at the west boundary.
    yoff: float, optional
        Latitude label offset from the axis (units are points).
        If negative (by default), the labels are drawn at the south boundary,
        otherwise at the north boundary.
    degree: bool, optional
        Add a degree symbol to tick labels.
    **text_kw: dict, optional
        Label text properties.
    """
    # Define what boundary to use
    extent = ax.get_extent(crs=ccrs.PlateCarree())
    if xoff <= 0:
        xpos = extent[0]  # labels at the east boundary
    else:
        xpos = extent[1]  # labels at the west boundary
    if yoff <= 0:
        ypos = extent[2]  # labels at the south boundary
    else:
        ypos = extent[3]  # labels at the north boundary
    geodetic_trans = ccrs.Geodetic()
    xlab_kw = ylab_kw = {**{"va": "center", "ha": "center"}, **text_kw}
    for xtick in xticks:
        s = fmt_lonlat(xtick, "lon", degree=degree)
        text_transform = offset_copy(
            geodetic_trans._as_mpl_transform(ax), fig=fig, units="points", x=0, y=yoff
        )
        ax.text(xtick, ypos, s, transform=text_transform, **xlab_kw)

    for ytick in yticks:
        s = fmt_lonlat(ytick, "lat", degree=degree)
        text_transform = offset_copy(
            geodetic_trans._as_mpl_transform(ax), fig=fig, units="points", x=xoff, y=0
        )
        ax.text(xpos, ytick, s, transform=text_transform, **ylab_kw)
