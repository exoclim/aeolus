# -*- coding: utf-8 -*-
"""Test coord submodule."""
from pathlib import Path

from aeolus import coord

import iris.coords
import iris.cube
import iris.exceptions

import numpy as np
import numpy.testing as npt

import pytest


iris.FUTURE.datum_support = True
TST_DATA = Path(__file__).parent / "data" / "test_data"


@pytest.fixture(scope="module")
def example_cubelist_from_file():
    return iris.load(str(TST_DATA / "netcdf" / "held_suarez_1200d.nc"))


def test_cell_centres():
    arr = np.array([-1, 0, 1, 2])
    des = np.array([-0.5, 0.5, 1.5])
    act = coord._cell_centres(arr)
    npt.assert_allclose(act, des)
    des = np.array([-0.8, 0.2, 1.2])
    act = coord._cell_centres(arr, bound_position=0.2)
    npt.assert_allclose(act, des)


def test_cell_bounds():
    arr = np.array([26.0, 27.0, 28.0, 29.0, 30.0])
    des = np.array([25.5, 26.5, 27.5, 28.5, 29.5, 30.5])
    act = coord._cell_bounds(arr)
    npt.assert_allclose(act, des)
    des = np.array([25.0, 27.0, 28.0, 29.0, 30.0, 31.0])
    act = coord._cell_bounds(arr, bound_position=1)
    npt.assert_allclose(act, des)


def test__is_longitude_global():
    assert coord._is_longitude_global(np.arange(0, 360, 2.5))
    assert coord._is_longitude_global(np.linspace(-180, 180, 10))
    assert not coord._is_longitude_global(np.arange(0, 180, 2.5))


def test_coord_to_cube():
    xc = iris.coords.DimCoord([-1, 2, 3], units="m", standard_name="longitude")
    yc = iris.coords.DimCoord([10, 30, 50, 70], units="m", standard_name="latitude")
    zc = iris.coords.DimCoord([1000, 500], units="hPa", standard_name="air_pressure")
    arr = np.arange(24).reshape((2, 4, 3))
    cube = iris.cube.Cube(
        data=arr,
        dim_coords_and_dims=[i[::-1] for i in [*enumerate((zc, yc, xc))]],
        standard_name="x_wind",
        units="m/s",
    )
    cc_bc = coord.coord_to_cube(cube, "longitude")
    assert cc_bc.ndim == 3
    assert cc_bc.coord("latitude") == cube.coord("latitude")
    cc = coord.coord_to_cube(cube, "longitude", broadcast=False)
    assert cc.ndim == 1
    assert cc.shape == xc.shape
    assert cc.standard_name == "longitude"
    assert cc.units == "m"
    npt.assert_allclose(cc.data, xc.points)


def test_isel_cube(example_cubelist_from_file):
    cube = example_cubelist_from_file[0]
    result = coord.isel(cube, "latitude", 11)
    assert result.shape == (7, 72)
    assert result.coord("latitude").points[0] == -44.0
    kgo = [
        279.34421964,
        304.32601524,
        314.35886024,
        357.66854755,
        443.61117028,
        575.30497115,
        732.65202309,
    ]
    npt.assert_allclose(result.data[:, 0], kgo)
    result = coord.isel(cube, cube.coords()[-1], 0)
    assert result.shape == (45, 72)
    with pytest.raises(IndexError):
        coord.isel(cube, "latitude", 123, skip_not_found=True)
    with pytest.raises(iris.exceptions.CoordinateNotFoundError):
        coord.isel(cube, "blah", 11)
    result = coord.isel(cube, "blah", 11, skip_not_found=True)
    assert result is cube


def test_isel_cubelist(example_cubelist_from_file):
    cl = example_cubelist_from_file
    result = coord.isel(cl, "latitude", 11)
    assert all([i.shape == (7, 72) for i in result])
    assert result[-1].coord("latitude").points[0] == -44.0
    result = coord.isel(cl, "blah", 11)
    assert result == cl
    with pytest.raises(iris.exceptions.CoordinateNotFoundError):
        result = coord.isel(cl, "blah", 11, skip_not_found=False)
    cube1 = cl[0].copy()
    cube1.remove_coord("latitude")
    cube1.remove_coord("longitude")
    cube2 = cl[1]
    result = coord.isel([cube1, cube2], "latitude", 11)
    assert len(result) == 2
