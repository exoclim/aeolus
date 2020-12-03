"""Metadata-related functionality."""
from collections.abc import Iterable
import functools

import cf_units

import iris


def update_metadata(name=None, units=None):
    """Update metadata of a cube returned by a function."""

    def update_name_units(cube):
        """Update name and convert units."""
        if isinstance(name, str):
            cube.rename(name)
        if isinstance(units, (str, cf_units.Unit)):
            cube.convert_units(units)

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Call the decorated function
            out = func(*args, **kwargs)
            if isinstance(out, iris.cube.Cube):
                update_name_units(out)
            elif isinstance(out, Iterable):
                [update_name_units(cube) for cube in out]
            return out

        return wrapper

    return decorator
