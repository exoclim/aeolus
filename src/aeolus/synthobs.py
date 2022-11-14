# -*- coding: utf-8 -*-
"""Functions for calculating synthetic observations."""
import warnings

from iris.analysis import SUM
from iris.coords import AuxCoord, DimCoord
from iris.cube import Cube
from iris.exceptions import CoordinateNotFoundError as CoNotFound

import numpy as np

from .coord import roll_cube_pm180
from .model import um


__all__ = (
    "calc_geom_mean_mirrored",
    "calc_stellar_flux",
    "calc_transmission_spectrum",
    "calc_transmission_spectrum_day_night_average",
    "read_normalized_stellar_flux",
    "read_spectral_bands",
)


def calc_geom_mean_mirrored(cube_a, cube_b, model=um):
    """
    Calculate geometric mean of two cubes with one of them rolled along the x-axis.

    This function can be used to get an average transmission flux
    calculated separately from the day- and night-side perspective.

    cube_a: iris.cube.Cube
        Cube with an x-coordinate.
    cube_b: iris.cube.Cube
        Another cube with an x-coordinate to be flipped before being averaged with cube A.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.
    """
    # Get x-coordinate from the 1st cube.
    x_coord = cube_a.coord(model.x)

    # Roll the 2nd cube by 180 degrees
    # This is specific to the UM output
    cube_b_rolled = roll_cube_pm180(cube_b, model=model)
    cube_b_rolled.replace_coord(x_coord)

    # Calculate geometric mean of the two cubes
    geom_mean = (cube_a * cube_b_rolled) ** 0.5
    return geom_mean


def calc_stellar_flux(spectral_file, stellar_constant_at_1_au):
    """
    Calculate the stellar flux per spectral band.

    Parameters
    ----------
    spectral_file: pathlib.Path
        Path to the location of a SOCRATES spectral file.
    stellar_constant_at_1_au: float or iris.cube.Cube
        Stellar constant at 1 AU [W m-2].

    Returns
    -------
    iris.cube.Cube
        Stellar flux per spectral band [W m-2].
    """
    # Ensure that an input constant is an iris cube
    if not isinstance(stellar_constant_at_1_au, Cube):
        stellar_constant_at_1_au = Cube(
            stellar_constant_at_1_au,
            long_name="stellar_constant_at_1_au",
            units="W m-2",
        )
    normalized_stellar_flux = read_normalized_stellar_flux(spectral_file)
    stellar_flux = normalized_stellar_flux * stellar_constant_at_1_au
    stellar_flux.rename("stellar_flux")
    return stellar_flux


def calc_transmission_spectrum(
    trans_flux,
    spectral_file=None,
    stellar_constant_at_1_au=None,
    stellar_radius=None,
    planet_top_of_atmosphere=None,
    model=um,
):
    r"""
    Convert the model output of transmission flux to a planetary-stellar radius ratio.

    Parameters
    ----------
    trans_flux: iris.cube.Cube
        Transmission flux on spectral bands and optionally latitudes and longitudes.
        In the Met Office Unified Model this is STASH items 555, 556, 755, 756 in section 1.
    spectral_file: pathlib.Path
        Path to the location of a SOCRATES spectral file.
    stellar_constant_at_1_au: float or iris.cube.Cube
        Stellar constant at 1 AU [W m-2].
    stellar_radius: float or iris.cube.Cube
        Stellar radius [m].
    planet_top_of_atmosphere: float or iris.cube.Cube
        The extent of the planetary atmosphere [m].
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    iris.cube.Cube
        The ratio of the effective planetary radius to the stellar radius per spectral band [1].
        Spectral band centres [m] is attached as an auxiliary coordinate.

    Notes
    -------
    The transmission spectrum is the ratio of the effective planetary radius to the stellar
    radius calculated per spectral band:

    .. math::

        \frac{R_p (\nu)}{R_s} = \sqrt{(\frac{R_{p,TOA}}{R_s})^2 -
         \frac{\sum_{lat,lon}^{}F_{transmitted} (\nu)}{F_{stellar} (\nu)}}

    where
    :math:`R_p(\nu)` is the effective planetary radius,
    :math:`R_s` is the stellar radius,
    :math:`R_{p,TOA}` is the extent of the planetary atmosphere (which usually is the sum of
    the planetary radius and the height of the model domain),
    :math:`\sum_{lat,lon}^{}F_{transmitted}(\nu)` is the total transmitted flux,
    :math:`F_{stellar}(\nu)` is the stellar flux.
    """
    # Ensure that input constants are iris cubes
    if not isinstance(stellar_constant_at_1_au, Cube):
        stellar_constant_at_1_au = Cube(
            stellar_constant_at_1_au,
            long_name="stellar_constant_at_1_au",
            units="W m-2",
        )
    if not isinstance(stellar_radius, Cube):
        stellar_radius = Cube(
            stellar_radius,
            long_name="stellar_radius",
            units="m",
        )
    if not isinstance(planet_top_of_atmosphere, Cube):
        planet_top_of_atmosphere = Cube(
            planet_top_of_atmosphere,
            long_name="planet_top_of_atmosphere",
            units="m",
        )
    # Calculate stellar flux
    stellar_flux = calc_stellar_flux(spectral_file, stellar_constant_at_1_au)

    # Sum transmission flux over all latitudes and longitudes
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            trans_flux = trans_flux.collapsed([model.y, model.x], SUM)
    except CoNotFound:
        pass
    # trans_flux.rename("shortwave_transmission_flux")
    trans_flux.units = "W m-2"

    # Calculate the ratio of the total transmitted flux to the stellar flux
    for coord_dim, coord in enumerate(trans_flux.dim_coords):
        if coord.name().lower() in ["pseudo", "pseudo_level"]:
            coord.rename("spectral_band_index")
            break
    flux_ratio = trans_flux / stellar_flux

    # Calculate the ratio of the effective planetary radius to the stellar radius
    rp_eff_over_rs_squared = (planet_top_of_atmosphere / stellar_radius) ** 2 - flux_ratio
    # Make all negative values zero
    rp_eff_over_rs_squared = rp_eff_over_rs_squared.copy(
        data=rp_eff_over_rs_squared.core_data().clip(min=0.0)
    )
    rp_eff_over_rs = rp_eff_over_rs_squared ** (0.5)
    rp_eff_over_rs.rename("ratio_of_effective_planetary_radius_to_stellar_radius")

    # Find spectral band centers
    spectral_bands = read_spectral_bands(spectral_file)
    spectral_band_centres = 0.5 * (
        spectral_bands["lower_wavelength_limit"] + spectral_bands["upper_wavelength_limit"]
    )
    spectral_bands_coord = AuxCoord(
        spectral_band_centres,
        long_name="spectral_band_centres",
        units="m",
    )
    spectral_bands_ll = AuxCoord(
        spectral_bands["lower_wavelength_limit"],
        long_name="spectral_band_lower_limit",
        units="m",
    )
    spectral_bands_ul = AuxCoord(
        spectral_bands["upper_wavelength_limit"],
        long_name="spectral_band_upper_limit",
        units="m",
    )

    # Attach spectral bands to the resulting cube as an auxiliary coordinate
    rp_eff_over_rs.add_aux_coord(spectral_bands_coord, data_dims=(coord_dim,))
    rp_eff_over_rs.add_aux_coord(spectral_bands_ll, data_dims=(coord_dim,))
    rp_eff_over_rs.add_aux_coord(spectral_bands_ul, data_dims=(coord_dim,))

    return rp_eff_over_rs


def calc_transmission_spectrum_day_night_average(
    trans_flux_day,
    trans_flux_night,
    spectral_file=None,
    stellar_constant_at_1_au=None,
    stellar_radius=None,
    planet_top_of_atmosphere=None,
    model=um,
):
    r"""
    Convert the model output of transmission flux to a planetary-stellar radius ratio.

    For UM output, this function averages the flux calculated from the day-side and the night-side
    of the planet. Why does it use a geometric mean? The reason to use a geometric average instead
    of an arithmetic average is that you want to add the optical depths. The flux decreases via
    Beer's law (i.e., it's proportional to :math:`exp[-optical depth]`) so when you multiply the
    dayside fluxes and nightside fluxes together and take a square root, you end up with
    :math:`exp[-( optical depth 1 + optical depth 2)/2]`. Since each optical depth is double the
    optical depth for it's respective side, the factors of two cancel and you end up with
    :math:`exp[-(true optical depth)]`.

    Parameters
    ----------
    trans_flux_day: iris.cube.Cube
        Transmission flux on spectral bands and optionally latitudes and longitudes.
        Day-side perspective.
        In the Met Office Unified Model this is STASH items 555, 556, 755, 756 in section 1.
    trans_flux_night: iris.cube.Cube
        Transmission flux on spectral bands and optionally latitudes and longitudes.
        Night-side perspective.
        In the Met Office Unified Model this is STASH items 555, 556, 755, 756 in section 1.

    For other parameters, see the docstring of `aeolus.synthobs.calc_transmission_spectrum()`

    See Also
    --------
    aeolus.synthobs.calc_transmission_spectrum
    """
    # Average the day and night flux
    trans_flux = calc_geom_mean_mirrored(trans_flux_day, trans_flux_night, model=model)
    return calc_transmission_spectrum(
        trans_flux,
        spectral_file=spectral_file,
        stellar_constant_at_1_au=stellar_constant_at_1_au,
        stellar_radius=stellar_radius,
        planet_top_of_atmosphere=planet_top_of_atmosphere,
        model=um,
    )


def read_normalized_stellar_flux(spectral_file):
    """
    Read the normalized stellar flux per spectral band from a SOCRATES spectral file.

    Parameters
    ----------
    spectral_file: pathlib.Path
        Path to the location of the SOCRATES spectral file.

    Returns
    -------
    iris.cube.Cube
        Normalized stellar flux per spectral band [1].
    """
    with spectral_file.open("r") as f:
        lines = []
        band_block = False
        for line in f:
            if line.startswith("*BLOCK: TYPE =    2"):
                band_block = True
            if band_block:
                if line.startswith("*END"):
                    break
                elif line.lstrip().split(" ")[0].isnumeric():
                    lines.append(
                        tuple([float(i) for i in line.lstrip().rstrip("\n").split(" ") if i])
                    )
    normalized_stellar_flux_arr = np.array(
        lines, dtype=[("spectral_band_index", "u4"), ("normalized_stellar_flux", "f4")]
    )
    # Compose an iris cube
    spectral_band_index = DimCoord(
        normalized_stellar_flux_arr["spectral_band_index"],
        long_name="spectral_band_index",
        units="1",
    )
    normalized_stellar_flux = Cube(
        normalized_stellar_flux_arr["normalized_stellar_flux"],
        long_name="normalized_stellar_flux",
        dim_coords_and_dims=[(spectral_band_index, 0)],
        units="1",
    )
    return normalized_stellar_flux


def read_spectral_bands(spectral_file):
    """
    Read spectral bands from a SOCRATES spectral file.

    Parameters
    ----------
    spectral_file: pathlib.Path
        Path to the location of a SOCRATES spectral file.

    Returns
    -------
    numpy.ndarray
        An array with a list of tuples describing spectral bands.
        Tuple structure:
        (spectral_band_index, lower_wavelength_limit, upper_wavelength_limit).
        Units [m] as stated in a spectral file but not checked directly.
    """
    with spectral_file.open("r") as f:
        lines = []
        band_block = False
        for line in f:
            if line.startswith("*BLOCK: TYPE =    1"):
                band_block = True
            if band_block:
                if line.startswith("*END"):
                    break
                elif line.lstrip().split(" ")[0].isnumeric():
                    lines.append(
                        tuple([float(i) for i in line.lstrip().rstrip("\n").split(" ") if i])
                    )
    spectral_bands = np.array(
        lines,
        dtype=[
            ("spectral_band_index", "u4"),
            ("lower_wavelength_limit", "f4"),
            ("upper_wavelength_limit", "f4"),
        ],
    )
    return spectral_bands
