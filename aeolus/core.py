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

    def __init__(
        self,
        files=None,
        name="",
        description="",
        planet="",
        const_dir=None,
        model_type=None,
        parent=None,
        children=None,
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
        parent: aeolus.core.Run, optional
            Pointer to this run's driving model if this is a LAM-type simulation.
        children: list, optional
            List of `aeolus.core.Run` objects if this is a driving model.

        See also
        --------
        aeolus.const.init_const
        """
        self.name = name
        self.description = description

        # Planetary constants
        self.const = init_const(planet, directory=const_dir)
        try:
            self.const.radius.convert_units("m")
            self._coord_system = iris.coord_systems.GeogCS(semi_major_axis=self.const.radius.data)
        except AttributeError:
            self._coord_system = None
            warn("Run initialised without a coordinate system.", AeolusWarning)

        # If the model is global or LAM (nested) and what its driving model is
        self.model_type = model_type
        self.parent = parent
        self.children = children

        if files is not None:
            self.load_data(files)
            try:
                cube_yx = self.raw.extract(DIM_CONSTR_YX_R)[0]
                self.domain = Region.from_cube(cube_yx, name=f"{name}_domain", shift_lons=True)
            except IndexError:
                warn("Run initialised without a domain.")

    def load_data(self, files):
        """Load cubes."""
        if isinstance(files, (list, set, tuple)):
            fnames = [str(i) for i in files]
        elif isinstance(files, (str, Path)):
            fnames = str(files)
        else:
            raise ArgumentError(f"Input type {type(files)} is not allowed.")
        self.raw = iris.load(fnames)

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
        self.proc = iris.cube.CubeList()
        if callable(func):
            self.proc = func(self.raw, **func_args)
        for cube in self.proc:
            # add constants to cube attributes
            cube.attributes["planet_conf"] = self.const
            for coord in cube.coords():
                if coord.coord_system:
                    # Replace coordinate system with the planet radius given in `self.const`
                    coord.coord_system = self._coord_system

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
