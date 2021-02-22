"""Core submodule of aeolus package."""
from warnings import warn

import iris
from iris.exceptions import ConstraintMismatchError as ConMisErr

from .calc import diag
from .calc.meta import copy_doc
from .const import add_planet_conf_to_cubes, init_const
from .coord import CoordContainer
from .exceptions import AeolusWarning
from .io import load_data, save_cubelist
from .model import um
from .region import Region
from .subset import DimConstr

__all__ = (
    "AtmosFlow",
    "Run",
)


class AtmosFlow:
    """
    Atmospheric Flow.

    Used to store and calculate atmospheric fields from gridded model output.
    Derived quantities are stored as cached properties to save computational
    time.

    Uses spherical coordinates.

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
        processed=False,
    ):
        """
        Instantiate an `AtmosFlow` object.

        Parameters
        ----------
        cubes: iris.cube.CubeList
            Atmospheric fields.
        name: str, optional
            The name of this `AtmosFlow`.
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

        See also
        --------
        aeolus.const.init_const, aeolus.core.Run
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
            warn("Initialised without a domain.", AeolusWarning)

        # Common coordinates
        self.coord = CoordContainer(self._cubes)

        # Variables as attributes
        self.assign_fields()

    def __getitem__(self, key):
        """Redirect self[key] to self.key."""
        return self.__getattribute__(key)

    def _assign_fields(self):
        """Assign input cubelist items as attributes of this class."""
        kwargs = {}
        for key in self.model.__dataclass_fields__:
            try:
                kwargs[key] = self._cubes.extract_cube(getattr(self.model, key))
            except ConMisErr:
                pass
        self.__dict__.update(**kwargs)
        del kwargs, key

    def _update_planet(self, planet="", const_dir=None):
        """Add or update planetary constants."""
        self.planet = planet
        self.const = init_const(planet, directory=const_dir)
        try:
            self.const.radius.convert_units("m")
            self._coord_system = iris.coord_systems.GeogCS(semi_major_axis=self.const.radius.data)
        except AttributeError:
            self._coord_system = None
            warn("Run initialised without a coordinate system.", AeolusWarning)

    def _add_planet_conf_to_cubes(self):
        """Add or update planetary constants container to cube attributes."""
        add_planet_conf_to_cubes(self._cubes, self.const)

    @property
    @copy_doc(diag.wind_speed)
    def wind_speed(self):
        cmpnts = []
        for cmpnt_key in ["u", "v", "w"]:
            try:
                cmpnts.append(self[cmpnt_key])
            except AttributeError:
                pass
        return diag.wind_speed(*cmpnts)


class Run:
    """
    A single model 'run', i.e. simulation.

    Attributes
    ----------
    name: str
        The run's name.
    const: aeolus.const.ConstContainer
        Physical constants used in calculations for this run.
    model: aeolus.model.Model, optional
        Model class with relevant coordinate names.
    """

    attr_keys = ["name", "planet", "model_type", "timestep"]

    def __init__(
        self,
        files=None,
        name="",
        planet="",
        const_dir=None,
        model=um,
        model_type=None,
        timestep=None,
        processed=False,
    ):
        """
        Instantiate a `Run` object.

        Parameters
        ----------
        files: str or pathlib.Path, optional
            Wildcard for loading files.
        name: str, optional
            The run's name.
        planet: str, optional
            Planet configuration. This is used to get appropriate physical constants.
            If not given, Earth physical constants are initialised.
        const_dir: pathlib.Path, optional
            Path to a folder with JSON files containing constants for a specific planet.
        model: aeolus.model.Model, optional
            Model class with relevant coordinate and variable names.
        model_type: str, optional
            Type of the model run, global or LAM.
        timestep: int, optional
            Model time step in s.
        parent: aeolus.core.Run, optional
            Pointer to this run's driving model if this is a LAM-type simulation.
        children: list, optional
            List of `aeolus.core.Run` objects if this is a driving model.
        processed: bool, optional
            If True, data from `files` is assigned to `proc` attribute.

        See also
        --------
        aeolus.const.init_const
        """
        warn(
            "Run is deprecated and will be removed in the next release. "
            "Use iris.cube.CubeList or AtmosFlow instead.",
            AeolusWarning,
        )
        self.name = name

        # Planetary constants
        self._update_planet(planet=planet, const_dir=const_dir)

        # Model-specific names of variables and coordinates
        self.model = model
        self.dim_constr = DimConstr(model=self.model)

        # If the model is global or LAM (nested) and what its driving model is
        self.model_type = model_type
        self.timestep = timestep
        self.processed = processed

        if files:
            self.load_data(files)
            try:
                if self.processed:
                    cube_yx = self.proc.extract(self.dim_constr.relax.yx)[0]
                else:
                    cube_yx = self.raw.extract(self.dim_constr.relax.yx)[0]
                self.domain = Region.from_cube(cube_yx, name=f"{name}_domain", shift_lons=True)
            except IndexError:
                warn("Run initialised without a domain.")
        else:
            warn("Run initialised without input files")

    def load_data(self, files):
        """Load cubes."""
        if self.processed:
            self.proc = load_data(files)
            self._add_planet_conf_to_cubes()
        else:
            self.raw = load_data(files)

    def _update_planet(self, planet="", const_dir=None):
        """Add or update planetary constants."""
        self.planet = planet
        self.const = init_const(planet, directory=const_dir)
        try:
            self.const.radius.convert_units("m")
            self._coord_system = iris.coord_systems.GeogCS(semi_major_axis=self.const.radius.data)
        except AttributeError:
            self._coord_system = None
            warn("Run initialised without a coordinate system.", AeolusWarning)

    def _add_planet_conf_to_cubes(self):
        """Add or update planetary constants container to cube attributes."""
        add_planet_conf_to_cubes(self._cubes, self.const)

    def proc_data(self, func=None, **func_args):
        """
        Post-process data for easier analysis and store it in `self.proc` attribute.

        Parameters
        ----------
        func: callable
            Function that takes `iris.cube.CubeList` as its first argument.
        **func_args: dict-like, optional
            Keyword arguments passed to `func`.
        """
        if self.processed:
            warn("Run's data is already processed. Skipping.", AeolusWarning)
        else:
            self.proc = iris.cube.CubeList()
            if callable(func):
                self.proc = func(self.raw, **func_args)
            self._add_planet_conf_to_cubes()
            self.processed = True

    def add_data(self, func=None, **func_args):
        """
        Calculate additional diagnostics (of type `iris.cube.Cube`) and add them to `self.proc`.

        Parameters
        ----------
        func: callable
            Function that takes `iris.cube.CubeList` (`self.proc`) as its first argument
            and appends new cubes to it (and does not return anything).
        **func_args: dict-like, optional
            Keyword arguments passed to `func`.
        """
        if callable(func):
            func(self.proc, **func_args)

    def to_file(self, path):
        """
        Save `proc` cubelist to a file with appropriate metadata.

        Parameters
        ----------
        path: str or pathlib.Path
            File path.
        """
        run_attrs = {}
        for key in self.attr_keys:
            if getattr(self, key):
                run_attrs[key] = str(getattr(self, key))
        save_cubelist(self.proc, path, **run_attrs)
