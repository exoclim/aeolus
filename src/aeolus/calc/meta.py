"""Metadata-related functionality."""
import functools
from collections.abc import Iterable
from dataclasses import is_dataclass

import cf_units

import iris
from iris.analysis import _dimensional_metadata_comparison
from iris.util import broadcast_to_shape

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


def preserve_shape(func):
    """Preserve shape of the output array."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        cube_in = args[0]
        if isinstance(cube_in, iris.cube.Cube):
            orig_shape = cube_in.shape
        else:
            raise ArgumentError(
                "`preserve_shape` decorator requires the first argument"
                "of the function and its output to be a cube."
            )
        # Call the decorated function
        cube_out = func(*args, **kwargs)
        out_name = cube_out.name()
        cell_methods = cube_out.cell_methods
        dim_map = []
        for ndim, _ in enumerate(orig_shape):
            for pair in _dimensional_metadata_comparison(cube_out, cube_in)["not_equal"]:
                if ndim not in cube_in.coord_dims(pair[0]):
                    dim_map.append(ndim)
        bc_data = broadcast_to_shape(cube_out.data, orig_shape, sorted(set(dim_map)))
        cube_out = cube_in.copy(data=bc_data)
        cube_out.rename(out_name)
        cube_out.cell_methods = cell_methods
        return cube_out

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
