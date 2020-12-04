"""Metadata-related functionality."""
import functools
from collections.abc import Iterable

import cf_units

import iris


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


def copy_doc(original):
    """Copy docstring from another function."""

    def wrapper(func):
        func.__doc__ = original.__doc__
        return func

    return wrapper
