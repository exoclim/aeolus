# -*- coding: utf-8 -*-
"""Test io submodule."""
from pathlib import Path

from aeolus import io

import iris.coords
import iris.cube

import numpy as np
import numpy.testing as npt


TST_DATA = Path(__file__).parent / "data" / "test_data"


def test_create_dummy_cube():
    # using n_res
    cube_n72 = io.create_dummy_cube(n_res=72)
    assert cube_n72.shape == (108, 144)
    assert cube_n72.coord("longitude").units == "degrees"
    assert cube_n72.coord("longitude").circular
    npt.assert_allclose(cube_n72.coord("longitude").points[:3], [1.25, 3.75, 6.25])
    npt.assert_allclose(cube_n72.coord("latitude").points[:3], [-89.16666667, -87.5, -85.83333333])
    # using nlat and nlon
    cube = io.create_dummy_cube(nlat=108, nlon=144)
    assert cube_n72 == cube
    # using other grids
    assert io.create_dummy_cube(n_res=72, endgame=False, grid_type="b") == io.create_dummy_cube(
        n_res=72, endgame=True, grid_type="a"
    )


def test_load_vert_lev_theta():
    levs = io.load_vert_lev(TST_DATA / "vert_lev" / "vertlevs_L38_29t_9s_40km")
    assert isinstance(levs, np.ndarray)
    assert levs.shape == (39,)
    npt.assert_allclose(levs[0], 0)
    npt.assert_allclose(levs[-1], 39254.833576)
    levs2 = io.load_vert_lev(TST_DATA / "vert_lev" / "vertlevs_L38_29t_9s_40km", lev_type="theta")
    npt.assert_allclose(levs, levs2)


def test_load_vert_lev_rho():
    levs = io.load_vert_lev(TST_DATA / "vert_lev" / "vertlevs_L38_29t_9s_40km", lev_type="rho")
    assert levs.shape == (38,)
    npt.assert_allclose(levs[0], 9.9982061118072)
    npt.assert_allclose(levs[-1], 36081.76331548462)


def test_save_cubelist():
    xc = iris.coords.DimCoord([-1, 2, 3], units="m", standard_name="longitude")
    yc = iris.coords.DimCoord([10, 30, 50, 70], units="m", standard_name="latitude")
    zc = iris.coords.DimCoord([1000, 500], units="hPa", standard_name="air_pressure")
    arr = np.arange(24).reshape((2, 4, 3))
    cube = iris.cube.Cube(
        data=arr,
        dim_coords_and_dims=[i[::-1] for i in [*enumerate((zc, yc, xc))]],
        standard_name="x_wind",
        units="m/s",
        attributes={"title": "Test data"},
    )
    cl = iris.cube.CubeList([cube])
    fname = TST_DATA / "test_save_cubelist.nc"
    io.save_cubelist(cl, fname, dummy_attr="foo bar")
    # Test that only the original attributes are retained
    assert cl[0].attributes == {"title": "Test data"}
    # Test that the file exists
    assert fname.is_file()
    # Load and compare
    result = iris.load_cube(str(fname))
    result.units == "m/s"
    assert "title" in result.attributes
    assert "dummy_attr" in result.attributes
    npt.assert_allclose(result.data, cube.data)
