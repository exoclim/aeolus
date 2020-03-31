"""Core submodule of aeolus package."""
from pathlib import Path
from warnings import warn

import iris

from .const import init_const
from .exceptions import AeolusWarning, ArgumentError
from .region import Region
from .subset import DIM_CONSTR_YX_R

__all__ = ("Run",)


class Run:
    """
    A single model 'run', i.e. simulation.

    Attributes
    ----------
    name: str
        The run's name.
    description: str
        A description of the run.
    const: aeolus.const.ConstContainer
        Physical constants used in calculations for this run.
    """

    attr_keys = ["name", "description", "planet", "model_type", "timestep", "parent", "children"]

    def __init__(
        self,
        files=None,
        name="",
        description="",
        planet="",
        const_dir=None,
        model_type=None,
        timestep=None,
        parent=None,
        children=None,
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
        description: str, optional
            A description of the model. This is not used internally by
            aeolus; it is solely for the user's information.
        planet: str, optional
            Planet configuration. This is used to get appropriate physical constants.
            If not given, Earth physical constants are initialised.
        const_dir: pathlib.Path, optional
            Path to a folder with JSON files containing constants for a specific planet.
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
        self.name = name
        self.description = description

        # Planetary constants
        self._update_planet(planet=planet, const_dir=const_dir)

        # If the model is global or LAM (nested) and what its driving model is
        self.model_type = model_type
        self.timestep = timestep
        self.parent = parent
        self.children = children
        self.processed = processed

        if files:
            self.load_data(files)
            try:
                if self.processed:
                    cube_yx = self.proc.extract(DIM_CONSTR_YX_R)[0]
                else:
                    cube_yx = self.raw.extract(DIM_CONSTR_YX_R)[0]
                self.domain = Region.from_cube(cube_yx, name=f"{name}_domain", shift_lons=True)
            except IndexError:
                warn("Run initialised without a domain.")
        else:
            warn("Run initialised without input files")

    def load_data(self, files):
        """Load cubes."""
        if isinstance(files, (list, set, tuple)):
            fnames = [str(i) for i in files]
        elif isinstance(files, (str, Path)):
            fnames = str(files)
        else:
            raise ArgumentError(f"Input type {type(files)} is not allowed.")
        if self.processed:
            self.proc = iris.load(fnames)
            self._add_planet_conf_to_cubes()
        else:
            self.raw = iris.load(fnames)

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
        """Add or update planetary constants container to cube attributes within `self.proc`."""
        for cube in self.proc:
            # add constants to cube attributes
            cube.attributes["planet_conf"] = self.const
            for coord in cube.coords():
                if coord.coord_system:
                    # Replace coordinate system with the planet radius given in `self.const`
                    coord.coord_system = self._coord_system

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

    def to_netcdf(self, path):
        """
        Save `proc` cubelist to a netCDF file with appropriate metadata.

        Parameters
        ----------
        path: str or pathlib.Path
            File path.
        """
        run_attrs = {}
        for key in self.attr_keys:
            if getattr(self, key):
                run_attrs[key] = str(getattr(self, key))
        # Remove planet_conf attribute before saving
        out = iris.cube.CubeList()
        old_attrs = {}
        for cube in self.proc:
            old_attrs[cube.name()] = cube.attributes.copy()
            new_attrs = {**cube.attributes, **run_attrs}
            try:
                new_attrs.pop("planet_conf")
            except KeyError:
                pass
            cube.attributes = new_attrs
            out.append(cube)
        iris.save(out, str(path))
        # Restore original attributes
        for cube in self.proc:
            cube.attributes = old_attrs[cube.name()]
