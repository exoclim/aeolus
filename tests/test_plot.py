# -*- coding: utf-8 -*-
"""Test the plot submodule."""
from aeolus import plot

import iris.coords
import iris.cube

import pytest


@pytest.fixture(scope="module")
def example_cube_2d():
    xc = iris.coords.DimCoord([-1, 2, 3], units="degrees", standard_name="longitude")
    yc = iris.coords.DimCoord([10, 30, 50, 70], units="degrees", standard_name="latitude")
    arr = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]
    cube = iris.cube.Cube(
        data=arr,
        dim_coords_and_dims=[i[::-1] for i in [*enumerate((yc, xc))]],
        standard_name="x_wind",
        units="m/s",
    )
    return cube


def test_cube_minmeanmax_str(example_cube_2d):
    assert plot.cube_minmeanmax_str(example_cube_2d) == "min=0 | mean=4 | max=11"


def test_unit_format():
    assert plot.unit_format(-1.234e-5) == "$-1.2\\times10^{-5}$"
    assert plot.unit_format(6) == "$6.0$"
    assert plot.unit_format(1.234e5, unit="m s^{-1}") == "$1.2\\times10^{5}$ $m$ $s^{-1}$"
    assert plot.unit_format(1.234, decimal_digits=1, precision=3) == "$1.200$"
    assert plot.unit_format(1.234, decimal_digits=3, precision=1) == "$1.2$"
