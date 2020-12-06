"""Metadata-related functionality."""
import functools
from collections.abc import Iterable
from dataclasses import is_dataclass

import cf_units

import iris

from ..exceptions import ArgumentError


def const_from_attrs(func):
    """Get constants container from the input cube attributes if not passed explicitly."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        const = kwargs.pop("const", None)
        if const is None:
            for arg in args:
                if isinstance(arg, iris.cube.Cube):
                    const = arg.attributes.get("planet_conf")
                elif isinstance(arg, iris.cube.CubeList):
                    try:
                        const = arg[0].attributes.get("planet_conf")
                    except IndexError:
                        const = None
        if is_dataclass(const):
            kwargs.update(const=const)
        else:
            raise ArgumentError(
                "Constants dataclass has to be an argument or in the cube attributes"
            )
        # Call the decorated function
        out = func(*args, **kwargs)
        return out

    return wrapper


def copy_doc(original):
    """Copy docstring from another function."""

    def wrapper(func):
        func.__doc__ = original.__doc__
        return func

    return wrapper


def update_metadata(name=None, units=None, attrs=None):
    """Update metadata of a cube returned by a function."""

    def _update(cube):
        """Update name, convert units and update attributes."""
        if isinstance(name, str):
            cube.rename(name)
        if isinstance(units, (str, cf_units.Unit)):
            cube.convert_units(units)
        if isinstance(attrs, dict):
            cube.attributes.update(attrs)

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Call the decorated function
            out = func(*args, **kwargs)
            if isinstance(out, iris.cube.Cube):
                _update(out)
            elif isinstance(out, Iterable):
                [_update(cube) for cube in out]
            return out

        return wrapper

    return decorator
