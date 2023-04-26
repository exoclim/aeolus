# -*- coding: utf-8 -*-
"""Contains functionality for making matplotlib plots."""
from .cart import GeoAxesGrid, label_global_map_gridlines
from .cloud import CloudPlotter
from .cm_custom import cloudtypes_cmap
from .mpl import (
    add_custom_legend,
    capitalise,
    figsave,
    hcross,
    linspace_pm1,
    make_list_2d,
    map_scatter,
    timeseries_1d,
    timeseries_2d,
)
from .text import (
    cube_minmeanmax_str,
    fmt_lonlat,
    subplot_label_generator,
    tex2cf_units,
    unit_format,
)


__all__ = (
    "add_custom_legend",
    "capitalise",
    "CloudPlotter",
    "cloudtypes_cmap",
    "cube_minmeanmax_str",
    "figsave",
    "fmt_lonlat",
    "GeoAxesGrid",
    "hcross",
    "label_global_map_gridlines",
    "linspace_pm1",
    "make_list_2d",
    "map_scatter",
    "subplot_label_generator",
    "tex2cf_units",
    "timeseries_1d",
    "timeseries_2d",
    "unit_format",
)
