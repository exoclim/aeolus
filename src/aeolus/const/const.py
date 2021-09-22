# -*- coding: utf-8 -*-
"""Main interface to the physical constants store."""
import json
from dataclasses import make_dataclass
from pathlib import Path

import iris.fileformats
from iris.coord_systems import GeogCS
from iris.cube import Cube

import numpy as np

from ..exceptions import ArgumentError, LoadError, _warn


__all__ = ("add_planet_conf_to_cubes", "get_planet_radius", "init_const")

CONST_DIR = Path(__file__).parent / "store"

DERIVED_CONST = {
    "dry_air_gas_constant": (
        lambda slf: slf.molar_gas_constant / slf.dry_air_molecular_weight,
        "J kg-1 K-1",
    ),
    "molecular_weight_ratio": (
        lambda slf: slf.condensible_molecular_weight / slf.dry_air_molecular_weight,
        "1",
    ),
    "poisson_exponent": (
        lambda slf: slf.dry_air_gas_constant / slf.dry_air_spec_heat_press,
        "1",
    ),
    "planet_rotation_rate": (lambda slf: (slf.day / (2 * np.pi)) ** (-1), "s-1"),
}


class ConstContainer:
    """Base class for creating dataclasses and storing planetary constants."""

    def __repr__(self):
        """Create custom repr."""
        cubes_str = ", ".join(
            [
                f"{getattr(self, _field).long_name} [{getattr(self, _field).units}]"
                for _field in self.__dataclass_fields__
            ]
        )
        return f"{self.__class__.__name__}({cubes_str})"

    def __post_init__(self):
        """Do things automatically after __init__()."""
        self._convert_to_iris_cubes()
        self._derive_const()

    def _convert_to_iris_cubes(self):
        """Loop through fields and convert each of them to `iris.cube.Cube`."""
        for name in self.__dataclass_fields__:
            _field = getattr(self, name)
            cube = Cube(data=_field.get("value"), units=_field.get("units", 1), long_name=name)
            object.__setattr__(self, name, cube)

    def _derive_const(self):
        """Not fully implemented yet."""
        for name, recipe in DERIVED_CONST.items():
            func, units = recipe
            try:
                cube = func(self)
                cube.convert_units(units)
                cube.rename(name)
                object.__setattr__(self, name, cube)
            except AttributeError:
                pass


def _read_const_file(name, directory=CONST_DIR):
    """Read constants from the JSON file."""
    if not isinstance(directory, Path):
        raise ArgumentError("directory must be a pathlib.Path object")
    try:
        with (directory / name).with_suffix(".json").open("r") as fp:
            list_of_dicts = json.load(fp)
        # transform the list of dictionaries into a dictionary
        const_dict = {}
        for vardict in list_of_dicts:
            const_dict[vardict["name"]] = {k: v for k, v in vardict.items() if k != "name"}
        return const_dict
    except FileNotFoundError:
        raise LoadError(
            f"JSON file for {name} configuration not found, check the directory: {directory}"
        )


def init_const(name="general", directory=None):
    """
    Create a dataclass with a given set of constants.

    Parameters
    ----------
    name: str, optional
        Name of the constants set.
        Should be identical to the JSON file name (without the .json extension).
        If not given, only general physical constants are returned.
    directory: pathlib.Path, optional
        Path to a folder with JSON files containing constants for a specific planet.

    Returns
    -------
    Dataclass with constants as iris cubes.

    Examples
    --------
    >>> c = init_const('earth')
    >>> c
    EarthConstants(gravity [m s-2], radius [m], day [s], solar_constant [W m-2], ...)
    >>> c.gravity
    <iris 'Cube' of gravity / (m s-2) (scalar cube)>
    """
    cls_name = f"{name.capitalize()}Constants"
    if directory is None:
        # use default directory
        kw = {}
    else:
        kw = {"directory": directory}
    # transform the list of dictionaries into a dictionary
    const_dict = _read_const_file("general")  # TODO: make this more flexible?
    if name != "general":
        const_dict.update(_read_const_file(name, **kw))
    kls = make_dataclass(
        cls_name,
        fields=[*const_dict.keys()],
        bases=(ConstContainer,),
        frozen=True,
        repr=False,
    )
    return kls(**const_dict)


def get_planet_radius(cube, default=iris.fileformats.pp.EARTH_RADIUS):
    """Get planet radius in metres from cube attributes or coordinate system."""
    cs = cube.coord_system("CoordSystem")
    if cs is not None:
        r = cs.semi_major_axis
    else:
        try:
            r = cube.attributes["planet_conf"].radius.copy()
            r.convert_units("m")
            r = float(r.data)
        except (KeyError, LoadError):
            _warn("Using default radius")
            r = default
    return r


def add_planet_conf_to_cubes(cubelist, const):
    """
    Add constants container to the cube attributes and replace its coordinate system.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        List of cubes containing a cube of zonal velocity (u).
    const: aeolus.const.const.ConstContainer, optional
        Constainer with the relevant planetary constants.
    """
    const.radius.convert_units("m")
    _coord_system = GeogCS(semi_major_axis=const.radius.data)
    for cube in cubelist:
        # add constants to cube attributes
        cube.attributes["planet_conf"] = const
        for coord in cube.coords():
            if coord.coord_system:
                # Replace coordinate system with the planet radius given in `self.const`
                coord.coord_system = _coord_system
