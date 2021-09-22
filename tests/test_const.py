# -*- coding: utf-8 -*-
"""Test calc submodule."""
from pathlib import Path

from aeolus import const
from aeolus.exceptions import ArgumentError, LoadError

import iris.cube

import numpy.testing as npt

import pytest


TST_DATA = Path(__file__).parent / "data" / "test_data"
CONST_FILE = "dummy"


def test_init_const_general():
    cnsts = const.init_const()
    assert str(cnsts).startswith("GeneralConstants")
    assert isinstance(cnsts, const.const.ConstContainer)
    for key in cnsts.__dataclass_fields__.keys():
        attr = getattr(cnsts, key)
        assert isinstance(attr, iris.cube.Cube)
        assert attr.ndim == 0  # XXX could be relaxed
    key = "stefan_boltzmann"
    assert key in cnsts.__dataclass_fields__
    cube = getattr(cnsts, key)
    assert cube.units == "W m-2 K-4"
    npt.assert_allclose(cube.data, 5.670367e-08)


def test_loaderror():
    with pytest.raises(LoadError):
        const.init_const(CONST_FILE)
    with pytest.raises(LoadError):
        const.init_const(CONST_FILE, directory=TST_DATA / "nonexistent_directory")


def test_argumenterror():
    with pytest.raises(ArgumentError):
        const.init_const(CONST_FILE, directory=str(TST_DATA / "json"))


def test_init_const_custom():
    cnsts = const.init_const(CONST_FILE, directory=TST_DATA / "json")
    assert str(cnsts).startswith("DummyConstants")
    assert isinstance(cnsts, const.const.ConstContainer)
    for key in cnsts.__dataclass_fields__.keys():
        attr = getattr(cnsts, key)
        assert isinstance(attr, iris.cube.Cube)
        assert attr.ndim == 0  # XXX could be relaxed
    key = "my_constant"
    assert key in cnsts.__dataclass_fields__
    cube = getattr(cnsts, key)
    assert cube.units == "m s-1"
    npt.assert_allclose(cube.data, 123)
