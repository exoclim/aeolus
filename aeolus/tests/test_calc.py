"""Test calc submodule."""
from aeolus import calc

from cf_units import Unit

import iris

import numpy as np
import numpy.testing as npt


def test_integrate():
    """Test integrate function."""
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
