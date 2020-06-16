# -*- coding: utf-8 -*-
"""Main interface to the physical constants store."""
import json
from dataclasses import make_dataclass
from pathlib import Path
from warnings import warn

import iris

import numpy as np

from ..exceptions import AeolusWarning, ArgumentError, LoadError


__all__ = ("init_const", "get_planet_radius")

CONST_DIR = Path(__file__).parent / "store"


class ScalarCube(iris.cube.Cube):
    """Cube without coordinates."""

    def __repr__(self):
        """Repr of this class."""
        return f"<ScalarCube of {self.long_name} [{self.units}]>"

    def __deepcopy__(self, memo):
        """Deep copy of a scalar cube."""
        return self.from_cube(self._deepcopy(memo))

    @property
    def asc(self):
        """Convert cube to AuxCoord for math ops."""
        return iris.coords.AuxCoord(
            np.asarray(self.data), units=self.units, long_name=self.long_name
        )

    @classmethod
    def from_cube(cls, cube):
        """Convert iris cube to ScalarCube."""
        return cls(**{k: getattr(cube, k) for k in ["data", "units", "long_name"]})


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
            cube = ScalarCube(
                data=_field.get("value"), units=_field.get("units", 1), long_name=name
            )
            object.__setattr__(self, name, cube)

    def _derive_const(self):
        """Not fully implemented yet."""
        derivatives = {
            "dry_air_gas_constant": lambda slf: slf.molar_gas_constant
            / slf.dry_air_molecular_weight,
            "molecular_weight_ratio": lambda slf: slf.condensible_molecular_weight
            / slf.dry_air_molecular_weight,
            "poisson_exponent": lambda slf: slf.dry_air_gas_constant / slf.dry_air_spec_heat_press,
        }
        for name, func in derivatives.items():
            try:
                cube = ScalarCube.from_cube(func(self))
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
        cls_name, fields=[*const_dict.keys()], bases=(ConstContainer,), frozen=True, repr=False
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
            warn("Using default radius", AeolusWarning)
            r = default
    return r
