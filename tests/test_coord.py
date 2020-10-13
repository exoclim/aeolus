"""Test coord submodule."""
from aeolus import coord

import numpy as np
import numpy.testing as npt


def test_cell_centres():
    """Test cell_centres."""
    arr = np.array([-1, 0, 1, 2])
    des = np.array([-0.5, 0.5, 1.5])
    act = coord._cell_centres(arr)
    npt.assert_allclose(act, des)
    des = np.array([-0.8, 0.2, 1.2])
    act = coord._cell_centres(arr, bound_position=0.2)
    npt.assert_allclose(act, des)


def test_cell_bounds():
    """Test cell_bounds."""
    arr = np.array([26.0, 27.0, 28.0, 29.0, 30.0])
    des = np.array([25.5, 26.5, 27.5, 28.5, 29.5, 30.5])
    act = coord._cell_bounds(arr)
    npt.assert_allclose(act, des)
    des = np.array([25.0, 27.0, 28.0, 29.0, 30.0, 31.0])
    act = coord._cell_bounds(arr, bound_position=1)
    npt.assert_allclose(act, des)


def test__is_longitude_global():
    """Test _is_longitude_global."""
    assert coord._is_longitude_global(np.arange(0, 360, 2.5))
    assert coord._is_longitude_global(np.linspace(-180, 180, 10))
    assert not coord._is_longitude_global(np.arange(0, 180, 2.5))
