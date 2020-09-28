"""Input and output functionality."""
import f90nml

import iris

import numpy as np


__all__ = ("load_multidir", "load_vert_lev")


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
    with path_to_file as nml_file:
        nml = f90nml.read(nml_file)
        levs = np.array(nml["vertlevs"][f"eta_{lev_type}"]) * nml["vertlevs"]["z_top_of_model"]
    return levs
