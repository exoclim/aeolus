# -*- coding: utf-8 -*-
"""Input and output functionality."""
from pathlib import Path

import iris
from iris.coords import AuxCoord
from iris.cube import Cube, CubeList
from iris.fileformats.pp import EARTH_RADIUS
from iris.fileformats.um import structured_um_loading
from iris.coord_systems import GeogCS

import numpy as np
from typing import Optional, Any, Sequence


__all__ = ("create_dummy_cube", "load_data", "load_multidir", "load_vert_lev", "save_cubelist")


def create_dummy_cube(
    nlat: Optional[int] = None,
    nlon: Optional[int] = None,
    n_res: Optional[int] = None,
    endgame: Optional[bool] = True,
    grid_type: Optional[str] = "a",
    pm180: Optional[bool] = False,
) -> Cube:
    """
    Create a dummy 2D cube with given resolution compatible with the UM grid.

    Parameters
    ----------
    nlat: int, optional
        Number of points in latitude. Not needed if `n_res` is given.
    nlon: int, optional
        Number of points in longitude. Not needed if `n_res` is given.
    n_res: int, optional
        N-notation resolution.
        If given, `nlat` and `nlon` are calculated automatically.
    endgame: bool, optional
        Use ENDGame grid.
        If True, "A" grid starts at lat=-90+dlat/2, lon=dlon/2.
        If False, "A" grid starts at lat=-90, lon=0.
        ENDGame's "A" grid is NewDyn "B" grid and vice versa.
    grid_type: str, optional
        Type of the UM grid.
          - A is main UM grid
          - B is staggered lorenz wind grid
          - Cu is staggered U wind grid
          - Cv is staggered V wind grid
    pm180: bool, optional
        Use -/+180 instead of 0 to 360 for the longitude span.

    Returns
    -------
    cube: iris.cube.Cube
        A dummy 2D cube of zeros.

    Examples
    --------
    >>> create_dummy_cube(nlat=90, nlon=144)
    >>> create_dummy_cube(n=96)
    """
    assert grid_type.lower() in ["a", "b", "cu", "cv"], f"Grid {grid_type} not valid."

    # Calculate the number of points given the N-notation resolution.
    if n_res is not None:
        nlon = 2 * n_res
        nlat = 3 * (n_res // 2)

    # Default grid set up is for ENDGame A Grid, equivalent to NewDyn B grid.
    if endgame:
        lat_shift = grid_type.lower() in ("b", "cv")
        lon_shift = grid_type.lower() in ("b", "cu")
    else:
        lat_shift = grid_type.lower() in ("a", "cu")
        lon_shift = grid_type.lower() in ("a", "cv")

    # Set coordinate system to match PP data.
    geog_cs = GeogCS(EARTH_RADIUS)

    # Latitude
    spacing = 180.0 / nlat
    if lat_shift:
        lower_bound = -90.0
        upper_bound = 90.0
        act_nlat = nlat + 1
    else:
        lower_bound = -90.0 + spacing / 2.0
        upper_bound = 90.0 - spacing / 2.0
        act_nlat = nlat
    lats = np.linspace(lower_bound, upper_bound, act_nlat)
    lat_coord = iris.coords.DimCoord(
        lats, standard_name="latitude", units="degrees", coord_system=geog_cs
    )
    lat_coord.guess_bounds()

    # Longitude
    spacing = 360.0 / nlon
    if lon_shift:
        lower_bound = 0.0
        upper_bound = 360.0 - spacing
    else:
        lower_bound = 0.0 + spacing / 2.0
        upper_bound = 360.0 - spacing / 2.0
    lons = np.linspace(lower_bound, upper_bound, nlon)
    if pm180:
        lons -= 180.0
    lon_coord = iris.coords.DimCoord(
        lons, standard_name="longitude", units="degrees", circular=True, coord_system=geog_cs
    )
    lon_coord.guess_bounds()

    # Create a numpy array of zeros.
    data = np.zeros((len(lats), len(lons)), dtype=np.int32)

    # Assemble the cube.
    cube = Cube(
        data,
        dim_coords_and_dims=((lat_coord, 0), (lon_coord, 1)),
        var_name="dummy_cube",
        units="1",
    )

    return cube


def load_data(files: Sequence, structured: Optional[bool] = False) -> CubeList:
    """Use `iris.load` with optional structured loading for PP files."""
    if structured:
        with structured_um_loading():
            cubes = iris.load(files)
    else:
        cubes = iris.load(files)
    return cubes


def load_multidir(
    path_mask: str, labels: Sequence[str], label_name: Optional[str] = "run"
) -> CubeList:
    """Load cubelists from multiple directories and merge."""
    joint_cl = CubeList()
    for label in labels:
        cl = iris.load(str(path_mask).format(label))
        for cube in cl:
            cube.attributes["um_version"] = ""  # FIXME
            cube.add_aux_coord(AuxCoord([label], long_name=label_name))
            joint_cl.append(cube)
    return joint_cl.merge()


def load_vert_lev(path_to_file: Path, lev_type: Optional[str] = "theta") -> np.array:
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


def save_cubelist(cubelist: CubeList, path: Path, **aux_attrs: Optional[Any]) -> None:
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
    out = CubeList()
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
