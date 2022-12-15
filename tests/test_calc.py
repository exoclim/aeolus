# -*- coding: utf-8 -*-
"""Test calc submodule."""
from pathlib import Path

from aeolus import calc
from aeolus.exceptions import AeolusWarning

from cf_units import Unit

import iris
import iris.coord_systems
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
    return iris.load(TST_DATA / "netcdf" / "held_suarez_1200d.nc")


def test_integrate():
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
    x_int = calc.integrate(cube, "longitude")
    npt.assert_allclose(x_int.data, np.array([[3.0, 15.0, 27.0, 39.0], [51.0, 63.0, 75.0, 87.0]]))
    assert x_int.units == Unit("m2/s")
    assert x_int.name() == "integral_of_x_wind_wrt_longitude"
    t_arr = np.array(
        [
            [-300000.0, -350000.0, -400000.0],
            [-450000.0, -500000.0, -550000.0],
            [-600000.0, -650000.0, -700000.0],
            [-750000.0, -800000.0, -850000.0],
        ]
    )
    p_int = calc.integrate(cube, "air_pressure")
    p_int.convert_units("W m-2")
    npt.assert_allclose(p_int.data, t_arr)


def test_spatial(example_cubelist_from_file):
    cube = example_cubelist_from_file[0]
    result = calc.spatial(cube, "max")
    assert isinstance(result, iris.cube.Cube)
    assert result.shape == (7,)
    assert len(result.coords()) == 8
    assert result.coord("longitude") == iris.coords.DimCoord(
        points=np.array([180.0], dtype=np.float32),
        bounds=np.array([[0.0, 360.0]]),
        standard_name="longitude",
        var_name="longitude",
        units="degrees",
        coord_system=iris.coord_systems.GeogCS(6371229.0),
    )
    assert result.cell_methods[0].method == "maximum"
    assert result.cell_methods[0].coord_names == ("latitude", "longitude")
    kgo = [
        312.88888771,
        315.46518042,
        321.68071322,
        395.68476945,
        490.07506524,
        613.5893357,
        754.39737563,
    ]
    npt.assert_allclose(result.data, kgo)
    result = calc.spatial(cube, "MIN")
    kgo = [
        259.68906506,
        273.29139143,
        298.464596,
        320.08245475,
        423.69745232,
        544.69610463,
        698.78197378,
    ]
    npt.assert_allclose(result.data, kgo)


def test_spatial_mean(example_cubelist_from_file):
    cube = example_cubelist_from_file[0]
    result = calc.spatial_mean(cube)
    assert result == calc.spatial(cube, "mean")
    kgo = [
        286.0241389900014,
        304.5367306861022,
        313.7161581899977,
        350.3582590713208,
        446.03368741766377,
        570.7474771484258,
        728.6475794834439,
    ]
    npt.assert_allclose(result.data, kgo)


def test_spatial_mean_twice(example_cubelist_from_file):
    cube = example_cubelist_from_file[0]
    result = calc.spatial_mean(cube)
    with pytest.warns(AeolusWarning):
        other = calc.spatial_mean(result)
    assert result == other


def test_time_mean_single():
    xc = iris.coords.DimCoord([0, 1], units="m", standard_name="longitude")
    tc = iris.coords.DimCoord(
        np.arange(10), units="hours since 1970-01-01 00:00:00", standard_name="time"
    )
    tc.guess_bounds()
    arr = np.arange(20).reshape((10, 2))
    cube = iris.cube.Cube(
        data=arr,
        dim_coords_and_dims=[i[::-1] for i in [*enumerate((tc, xc))]],
        standard_name="x_wind",
        units="m/s",
    )
    result = calc.time_mean(cube)
    assert result.shape == (2,)
    assert result.units == "m/s"
    assert result.standard_name == "x_wind"
    kgo = [9.0, 10.0]
    npt.assert_allclose(result.data, kgo)
    result2 = calc.time_mean(cube, squeeze=True)
    assert result2 == result


def test_time_mean_collection():
    xc = iris.coords.DimCoord([0, 1], units="m", standard_name="longitude")
    tc = iris.coords.DimCoord(
        np.arange(10), units="hours since 1970-01-01 00:00:00", standard_name="time"
    )
    tc.guess_bounds()
    arr = np.arange(20).reshape((10, 2))
    cube1 = iris.cube.Cube(
        data=arr,
        dim_coords_and_dims=[i[::-1] for i in [*enumerate((tc, xc))]],
        standard_name="x_wind",
        units="m/s",
    )
    cube2 = iris.cube.Cube(
        data=arr * 2,
        dim_coords_and_dims=[i[::-1] for i in [*enumerate((tc, xc))]],
        standard_name="air_pressure",
        units="Pa",
    )
    cont_cl = iris.cube.CubeList([cube1, cube2])
    result = calc.time_mean(cont_cl)
    assert isinstance(result, iris.cube.CubeList)
    assert result[0] == calc.time_mean(cube1)
    assert result[1] == calc.time_mean(cube2)
    cont_tuple = (cube1, cube2)
    result = calc.time_mean(cont_tuple)
    assert isinstance(result, tuple)
    assert result[0] == calc.time_mean(cube1)
    assert result[1] == calc.time_mean(cube2)
    cont_dict = {"a": cube1, "b": cube2}
    result = calc.time_mean(cont_dict)
    assert isinstance(result, dict)
    assert result["a"] == calc.time_mean(cube1)
    assert result["b"] == calc.time_mean(cube2)
