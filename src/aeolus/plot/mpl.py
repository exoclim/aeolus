# -*- coding: utf-8 -*-
"""Matplotlib-related utilities."""
from matplotlib.lines import Line2D

__all__ = ("add_custom_legend",)


def add_custom_legend(ax_or_fig, styles_and_labels, **leg_kw):
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
