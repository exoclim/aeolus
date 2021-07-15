# -*- coding: utf-8 -*-
"""Test the plot submodule."""
from aeolus import plot


def test_unit_format():
    assert plot.unit_format(-1.234e-5) == "$-1.2\\times10^{-5}$"
    assert plot.unit_format(6) == "$6.0$"
    assert plot.unit_format(1.234e5, unit="m s^{-1}") == "$1.2\\times10^{5}$ $m$ $s^{-1}$"
    assert plot.unit_format(1.234, decimal_digits=1, precision=3) == "$1.200$"
    assert plot.unit_format(1.234, decimal_digits=3, precision=1) == "$1.2$"
