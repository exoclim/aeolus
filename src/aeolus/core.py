# -*- coding: utf-8 -*-
"""Core submodule of aeolus package."""
from cached_property import cached_property

from iris.coord_systems import GeogCS
from iris.exceptions import ConstraintMismatchError as ConMisErr

from .calc import diag
from .const import add_planet_conf_to_cubes, init_const
from .coord import CoordContainer
from .decor import ReprAtmoSimBase
from .exceptions import _warn
from .meta import copy_doc
from .model import um
from .region import Region
from .subset import DimConstr

__all__ = ("AtmoSim", "AtmoSimBase")


class AtmoSimBase:
    """
    Base class for creating atmospheric model simulation classes in aeolus.

    Used to store and calculate atmospheric fields from gridded model output.
    Derived quantities are stored as cached properties to save computational time.

    Assumes the data are in spherical coordinates on a regular longitude-latitude grid.

    Attributes
    ----------
    name: str
        The run's name.
    description: str
        A description of the run.
    const: aeolus.const.ConstContainer
        Physical constants used in calculations for this run.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate and variable names.
    """

    def __init__(
        self,
        cubes=None,
        name="",
        description="",
        planet="",
        const_dir=None,
        model=um,
        model_type=None,
        timestep=None,
        vert_coord=None,
        coord_check=False,
    ):
        """
        Instantiate an `AtmoSimBase` object.

        Parameters
        ----------
        cubes: iris.cube.CubeList
            Atmospheric fields.
        name: str, optional
            The name or label of this `AtmoSim`.
        description: str, optional
            This is not used internally; it is solely for the user's information.
        planet: str, optional
            Planet configuration. This is used to get appropriate physical constants.
            If not given, Earth physical constants are initialised.
        const_dir: pathlib.Path, optional
            Path to a folder with files containing constants for a specific planet.
        model: aeolus.model.Model, optional
            Model class with relevant coordinate and variable names.
        model_type: str, optional
            Type of the model run, global or LAM.
        timestep: int, optional
            Model time step in s.
        vert_coord: str, optional
            Character identificator for the type of the model's vertical coordinate.
            "z" - data on "level_height"
            "p" - data on pressure levels
        coord_check: bool, optional
            Check if all cubes have the same set of coordinates.

        See also
        --------
        aeolus.const.init_const, aeolus.coord.CoordContainer
        """
        self._cubes = cubes
        self.name = name
        self.description = description

        # Planetary constants
        self._update_planet(planet=planet, const_dir=const_dir)
        self._add_planet_conf_to_cubes()  # TODO: remove?

        # Model-specific names of variables and coordinates
        self.model = model
        self.dim_constr = DimConstr(model=self.model)

        # If the model is global or LAM (nested) and what its driving model is
        self.model_type = model_type
        self.timestep = timestep

        # Domain
        try:
            cube_yx = self._cubes.extract(self.dim_constr.relax.yx)[0]
            self.domain = Region.from_cube(cube_yx, name=f"{name}_domain", shift_lons=True)
        except IndexError:
            _warn("Initialised without a domain.")
            self.domain = None

        # Variables as attributes
        self._assign_fields()

        # Common coordinates
        dim_seq = ["tyx", "yx"]
        self.vert_coord = vert_coord
        # TODO: make it more flexible
        if self.vert_coord == "z":
            _constr = self.dim_constr.relax.z
        elif self.vert_coord == "p":
            _constr = self.dim_constr.relax.p
        else:
            _constr = self.dim_constr.relax.yx
        self.coord = CoordContainer(self._cubes.extract(_constr), coord_check=coord_check)

        if self.vert_coord is not None:
            dim_seq += [f"t{self.vert_coord}yx", f"{self.vert_coord}yx"]
        for seq in dim_seq:
            try:
                setattr(
                    self,
                    f"_ref_{seq}",
                    self._cubes.extract(getattr(self.dim_constr.strict, seq))[0],
                )
            except IndexError:
                pass

        self.repr = ReprAtmoSimBase(self)

    def __repr__(self):  # noqa
        return self.repr.str_repr(short=True)

    def __str__(self):  # noqa
        return self.repr.str_repr(short=False)

    def _repr_html_(self):
        return self.repr.html_repr()

    @classmethod
    def from_parent_class(cls, obj):
        """Dynamically inherit from a similar class."""
        new_obj = cls(
            cubes=obj._cubes,
            name=obj.name,
            description=obj.description,
            planet=obj.planet,
            const_dir=obj.const_dir,
            model=obj.model,
            model_type=obj.model_type,
            timestep=obj.timestep,
            vert_coord=obj.vert_coord,
        )
        return new_obj

    def __getitem__(self, key):
        """Redirect self[key] to self.key."""
        return self.__getattribute__(key)

    def _assign_fields(self):
        """Assign input cubelist items as attributes of this class."""
        self.vars = []
        kwargs = {}
        for key in self.model.__dataclass_fields__:
            try:
                kwargs[key] = self._cubes.extract_cube(getattr(self.model, key))
                self.vars.append(key)
            except ConMisErr:
                pass
        self.__dict__.update(**kwargs)
        del kwargs, key

    def _update_planet(self, planet="", const_dir=None):
        """Add or update planetary constants."""
        self.planet = planet
        self.const_dir = const_dir
        self.const = init_const(self.planet, directory=self.const_dir)
        try:
            self.const.radius.convert_units("m")
            self._coord_system = GeogCS(semi_major_axis=self.const.radius.data)
        except AttributeError:
            self._coord_system = None
            _warn("Initialised without a coordinate system.")

    def _add_planet_conf_to_cubes(self):
        """Add or update planetary constants container to cube attributes."""
        add_planet_conf_to_cubes(self._cubes, self.const)

    def extract(self, constraints):
        """
        Subset `AtmoSim` using iris constraints.

        Parameters
        ----------
        constraints: iris.Constraint or iterable of constraints
            A single constraint or an iterable.
        """
        new_obj = self.__class__(
            cubes=self._cubes.extract(constraints),
            name=self.name,
            description=self.description,
            planet=self.planet,
            const_dir=self.const_dir,
            model=self.model,
            model_type=self.model_type,
            timestep=self.timestep,
            vert_coord=self.vert_coord,
        )
        return new_obj


class AtmoSim(AtmoSimBase):
    """
    Main class for dealing with a atmospheric model simulation output in aeolus.

    Used to store and calculate atmospheric fields from gridded model output.
    Derived quantities are stored as cached properties to save computational time.
    """

    @cached_property
    @copy_doc(diag.sigma_p)
    def sigma_p(self):
        # TODO: pass the cube only?
        return diag.sigma_p(self._cubes, model=self.model)

    @cached_property
    @copy_doc(diag.wind_speed)
    def wind_speed(self):
        cmpnts = []
        for cmpnt_key in ["u", "v", "w"]:
            try:
                cmpnts.append(self[cmpnt_key])
            except AttributeError:
                pass
        return diag.wind_speed(*cmpnts)

    @cached_property
    @copy_doc(diag.toa_net_energy)
    def toa_net_energy(self):
        return diag.toa_net_energy(self._cubes, model=self.model)

    @cached_property
    @copy_doc(diag.sfc_water_balance)
    def sfc_water_balance(self):
        return diag.sfc_water_balance(self._cubes, const=self.const, model=self.model)
