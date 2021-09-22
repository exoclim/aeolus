# -*- coding: utf-8 -*-
"""Test the meta submodule."""
from aeolus import meta
from aeolus.const import init_const
from aeolus.const.const import ConstContainer
from aeolus.exceptions import ArgumentError

import iris.cube

import pytest


def test_const_from_attrs_default():
    @meta.const_from_attrs()
    def foo(arg, const=None):
        """Test function."""
        return const

    constants = init_const()
    cube = iris.cube.Cube([], attributes={"planet_conf": constants})
    result = foo(cube)
    assert isinstance(result, ConstContainer)
    result = foo(iris.cube.CubeList([cube]))
    assert isinstance(result, ConstContainer)
    with pytest.raises(ArgumentError):
        foo(0)
    with pytest.raises(ArgumentError):
        foo(iris.cube.Cube([]))
    with pytest.raises(ArgumentError):
        foo(iris.cube.Cube([]), const=1234)
    with pytest.raises(ArgumentError):
        foo(iris.cube.Cube([], attributes={"foo": constants}))
    with pytest.raises(ArgumentError):
        foo(iris.cube.Cube([], attributes={"planet_conf": 1234}))


def test_const_from_attrs_relax():
    @meta.const_from_attrs(strict=False)
    def bar(arg, const=None):
        """Test function."""
        return const

    constants = init_const()
    cube = iris.cube.Cube([], attributes={"planet_conf": constants})
    result = bar(cube)
    assert isinstance(result, ConstContainer)
    result = bar(iris.cube.CubeList([cube]))
    assert isinstance(result, ConstContainer)
    # with pytest.warns(AeolusWarning):
    #     bar(iris.cube.Cube([]))
    # with pytest.warns(AeolusWarning):
    #     bar(iris.cube.Cube([]), const=1234)
    # with pytest.warns(AeolusWarning):
    #     bar(iris.cube.Cube([], attributes={"foo": constants}))
    # with pytest.warns(AeolusWarning):
    #     bar(iris.cube.Cube([], attributes={"planet_conf": 1234}))


def test_copy_doc():
    def foo():
        """Make docstring."""
        pass

    @meta.copy_doc(foo)
    def bar():
        pass

    assert foo.__doc__ == bar.__doc__
