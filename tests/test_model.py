# -*- coding: utf-8 -*-
"""Test model submodule."""

from aeolus.model import um, um_stash
from aeolus.model.base import Model


def test_model_um():
    assert isinstance(um, Model)
    assert isinstance(um.u, str)
    assert um.u == "x_wind"


def test_model_um_stash():
    assert isinstance(um_stash, Model)
    assert isinstance(um_stash.u, str)
    assert um_stash.u == "m01s00i002"


def test_model_dummy():
    dummy = Model(u="abc")
    assert isinstance(dummy.u, str)
    assert dummy.u == "abc"
    assert dummy.v is None
