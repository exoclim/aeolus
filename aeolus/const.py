# -*- coding: utf-8 -*-
"""Meteorological constants as collections of scalar `iris` cubes."""
import json
from dataclasses import make_dataclass
from pathlib import Path

import iris

from .exceptions import LoadError

CONST_DIR = Path(__file__).parent / "phys_const_store"


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

    def _convert_to_iris_cubes(self):
        """Loop through fields and convert each of them to `iris.cube.Cube`."""
        for name in self.__dataclass_fields__:
            _field = getattr(self, name)
            cube = iris.cube.Cube(
                data=_field.get("value"), units=_field.get("units", 1), long_name=name
            )
            object.__setattr__(self, name, cube)

    def _derive_const(self):
        """Not implemented."""
        pass


def _read_const_file(name, directory=CONST_DIR):
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
            f"JSON file for {directory} {name} configuration not found, check the directory"
        )


def init_const(name, directory=None):
    """
    Create a dataclass with a given set of constants.

    Parameters
    ----------
    name: str
        Name of the constants set. Should be identical to the JSON file name.
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
    const_dict.update(_read_const_file(name, **kw))
    kls = make_dataclass(
        cls_name, fields=[*const_dict.keys()], bases=(ConstContainer,), frozen=True, repr=False
    )
    return kls(**const_dict)
