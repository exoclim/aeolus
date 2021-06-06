"""Test coord submodule."""
from aeolus import coord

import iris

import numpy as np
import numpy.testing as npt


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
