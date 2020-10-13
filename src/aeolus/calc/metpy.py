"""Interface to metpy calc functions."""
import functools

import cf_units

import iris

import metpy.units as metunits

import numpy as np

import pint

import xarray as xr

from ..exceptions import UnitFormatError


__all__ = ("preprocess_iris",)


def preprocess_iris(f):
    """
    Wrap a function from `metpy.calc` for it to accept iris cubes as arguments.

    In addition, this decorator converts `metpy.calc` output by using the first input argument
    as a 'donor' cube.
    Note this works only for functions that preserve dimensions and may not work with some units.
    Now that metpy has xarray preprocessor, this decorator depends on it.
    """
    # Define local functions
    def to_xarray(cube):
        """Convert `iris.cube.Cube` to `xarray.DataArray` and format units correctly."""
        _unit = None
        for ut_format in set(cf_units.UT_FORMATS):
            try:
                _unit = metunits.units(cube.units.format(ut_format))
            except pint.errors.DimensionalityError:
                pass
        if _unit is None:
            raise UnitFormatError(f"Unable to convert cube units of\n{repr(cube)}\nto metpy units")
        arr = xr.DataArray.from_iris(cube)
        arr.attrs["units"] = str(_unit)
        return arr

    def to_iris(donor_cube, arr, name):
        """Convert metpy calc result to `iris.cube.Cube`."""
        try:
            data = arr.magnitude
            units = str(arr.units).replace(" ** ", "^").replace(" * ", " ")
        except AttributeError:
            # output is probably a numpy array without units
            data = np.asarray(arr)
            units = "1"
        cube_out = donor_cube.copy(data=data)
        try:
            cube_out.units = units
        except ValueError:
            cube_out.units = "unknown"
        cube_out.rename(name)
        cube_out.attributes.pop("STASH", None)
        return cube_out

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        nargs = []
        _cube = None  # a donor cube with metadata

        # Loop over cubes and re-format units
        for arg in args:
            if isinstance(arg, iris.cube.Cube):
                if arg.ndim > 0:
                    # TODO: make this flexible
                    _cube = arg
                elif _cube is None:
                    _cube = arg

                nargs.append(to_xarray(arg))
            else:
                nargs.append(arg)

        kwargs = {
            k: (to_xarray(v) if isinstance(v, iris.cube.Cube) else v) for k, v in kwargs.items()
        }

        # Call the decorated function
        out = f(*nargs, **kwargs)
        if _cube is None:
            # Inputs are not iris cubes
            return out

        # Convert ouput to iris
        # TODO: check masked array behaviour
        if isinstance(out, (tuple, list, set)):
            res = []
            for i, iout in enumerate(out):
                res.append(to_iris(_cube, iout, f"{f.__name__}_output_{i}"))
            res = tuple(res)
        else:
            res = to_iris(_cube, out, f.__name__)
        return res

    return wrapper
