"""Functions for calculating synthetic observations."""
import iris
from iris.util import reverse

import numpy as np

from .coord import roll_cube_pm180
from .model import um


__all__ = (
    "read_spectral_bands",
    "read_normalized_stellar_flux",
    "calc_stellar_flux",
    "calc_transmission_spectrum_day_night_average",
)


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
    spectral_band_index = iris.coords.DimCoord(
        normalized_stellar_flux_arr["spectral_band_index"],
        long_name="spectral_band_index",
        units="1",
    )
    normalized_stellar_flux = iris.cube.Cube(
        normalized_stellar_flux_arr["normalized_stellar_flux"],
        long_name="normalized_stellar_flux",
        dim_coords_and_dims=[(spectral_band_index, 0)],
        units="1",
    )
    return normalized_stellar_flux


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
    if not isinstance(stellar_constant_at_1_au, iris.cube.Cube):
        stellar_constant_at_1_au = iris.cube.Cube(
            stellar_constant_at_1_au,
            long_name="stellar_constant_at_1_au",
            units="W m-2",
        )
    normalized_stellar_flux = read_normalized_stellar_flux(spectral_file)
    stellar_flux = normalized_stellar_flux * stellar_constant_at_1_au
    stellar_flux.rename("stellar_flux")
    return stellar_flux


def calc_transmission_spectrum_day_night_average(
    spectral_file,
    stellar_constant_at_1_au,
    stellar_radius,
    planet_top_of_atmosphere,
    planet_transmission_day,
    planet_transmission_night,
    model=um,
):
    r"""
    Calculate the synthetic transmission spectrum averaged over the dayside and the nightside.

    Parameters
    ----------
    spectral_file: pathlib.Path
        Path to the location of a SOCRATES spectral file.
    stellar_constant_at_1_au: float or iris.cube.Cube
        Stellar constant at 1 AU [W m-2].
    stellar_radius: float or iris.cube.Cube
        Stellar radius [m].
    planet_top_of_atmosphere: float or iris.cube.Cube
        The extent of the planetary atmosphere [m].
    planet_transmission_day: iris.cube.Cube
        Met Office Unified Model output of STASH item m01s01i755
        (in the case of hot Jupiters) from the dayside calculation.
    planet_transmission_night: iris.cube.Cube
        Met Office Unified Model output of STASH item m01s01i755
        (in the case of hot Jupiters) from the nightside calculation.
    model: aeolus.model.Model, optional
        Model class with relevant variable names.

    Returns
    -------
    numpy.ndarray
        Spectral band centers [m].
    iris.cube.Cube
        The ratio of the effective planetary radius to the stellar radius per spectral band [1].

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
    if not isinstance(stellar_constant_at_1_au, iris.cube.Cube):
        stellar_constant_at_1_au = iris.cube.Cube(
            stellar_constant_at_1_au,
            long_name="stellar_constant_at_1_au",
            units="W m-2",
        )
    if not isinstance(stellar_radius, iris.cube.Cube):
        stellar_radius = iris.cube.Cube(
            stellar_radius,
            long_name="stellar_radius",
            units="m",
        )
    if not isinstance(planet_top_of_atmosphere, iris.cube.Cube):
        planet_top_of_atmosphere = iris.cube.Cube(
            planet_top_of_atmosphere,
            long_name="planet_top_of_atmosphere",
            units="m",
        )

    # Load UM output from the dayside calculation
    day = planet_transmission_day
    day_lon_coord = day.coord(um.x)

    # Load UM output from the nightside calculation
    # Roll nightside data by 180 degrees
    night_rolled = roll_cube_pm180(planet_transmission_night)
    # Reverse longitude order
    night = reverse(night_rolled, night_rolled.coord_dims(um.x))
    # Replace the longitude coordinate to be able to do maths with iris
    night.replace_coord(day_lon_coord)

    # Use the same name for the spectral band index coordinate
    for coord in day.coords():
        if coord.name() in ["pseudo", "pseudo_level"]:
            coord.rename("spectral_band_index")
    for coord in night.coords():
        if coord.name() in ["pseudo", "pseudo_level"]:
            coord.rename("spectral_band_index")

    # Calculate the geometric mean of the dayside and nightside transmitted flux
    # and sum this flux over all latitudes and longitudes
    transmitted_flux = ((day * night) ** (0.5)).collapsed([um.y, um.x], iris.analysis.SUM)
    transmitted_flux.rename("total_transmitted_flux")
    transmitted_flux.units = "W m-2"

    # Calculate stellar flux
    stellar_flux = calc_stellar_flux(spectral_file, stellar_constant_at_1_au)

    # Calculate the ratio of the total transmitted flux to the stellar flux
    flux_ratio = transmitted_flux.copy(data=transmitted_flux.core_data() / stellar_flux.core_data())
    flux_ratio.rename("flux_ratio")
    flux_ratio.units = transmitted_flux.units / stellar_flux.units

    # Calculate the ratio of the effective planetary radius to the stellar radius
    rp_eff_over_rs_squared = (planet_top_of_atmosphere / stellar_radius) ** 2 - flux_ratio
    rp_eff_over_rs_squared.data[rp_eff_over_rs_squared.data < 0.0] = 0.0
    rp_eff_over_rs = rp_eff_over_rs_squared ** (0.5)
    rp_eff_over_rs.rename("radius_ratio")

    # Find spectral band centers
    spectral_bands = read_spectral_bands(spectral_file)
    spectral_band_centers = 0.5 * (
        spectral_bands["lower_wavelength_limit"] + spectral_bands["upper_wavelength_limit"]
    )

    return spectral_band_centers, rp_eff_over_rs
