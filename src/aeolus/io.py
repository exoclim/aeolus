"""Input and output functionality."""
from pathlib import Path

import iris

import numpy as np

from .exceptions import ArgumentError


__all__ = ("load_data", "load_multidir", "load_vert_lev", "save_cubelist")


def load_data(files):
    """Wrap `iris.load` to deal with `pathlib.Path` objects."""
    if isinstance(files, (list, set, tuple)):
        fnames = [str(i) for i in files]
    elif isinstance(files, (str, Path)):
        fnames = str(files)
    else:
        raise ArgumentError(f"Input type {type(files)} is not allowed.")
    return iris.load(fnames)


def load_multidir(path_mask, labels, label_name="run"):
    """Load cubelists from multiple directories and merge."""
    joint_cl = iris.cube.CubeList()
    for label in labels:
        cl = iris.load(str(path_mask).format(label))
        for cube in cl:
            cube.attributes["um_version"] = ""  # FIXME
            cube.add_aux_coord(iris.coords.AuxCoord([label], long_name=label_name))
            joint_cl.append(cube)
    return joint_cl.merge()


def load_vert_lev(path_to_file, lev_type="theta"):
    """
    Read data from the UM vertical levels file.

    Parameters
    ----------
    path_to_file: pathlib.Path
        Full path to the vertical levels file.
    lev_type: str, optional
        What levels to return: "theta" or "rho".

    Returns
    -------
    levs: numpy.array
        Array of height levels.
    """
    import f90nml  # noqa

    with path_to_file.open("r") as nml_file:
        nml = f90nml.read(nml_file)
        levs = np.array(nml["vertlevs"][f"eta_{lev_type}"]) * nml["vertlevs"]["z_top_of_model"]
    return levs


def save_cubelist(cubelist, path, **aux_attrs):
    """
    Save a cubelist w/o the `planet_conf` container to a file.

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        Cube list to write to disk.
    path: str or pathlib.Path
        File path.
    aux_attrs: dict, optional
        Dictionary of additional attributes to save with the cubes.
    """
    # Remove planet_conf attribute before saving
    out = iris.cube.CubeList()
    old_attrs = []
    for cube in cubelist:
        old_attrs.append(cube.attributes.copy())
        new_attrs = {**cube.attributes, **aux_attrs}
        try:
            new_attrs.pop("planet_conf")
        except KeyError:
            pass
        cube.attributes = new_attrs
        out.append(cube)
    iris.save(out, str(path))
    # Restore original attributes
    for cube, attrs in zip(cubelist, old_attrs):
        cube.attributes = attrs
