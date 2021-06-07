"""Test synthobs submodule."""
from pathlib import Path

from aeolus import synthobs

import iris

import numpy as np
import numpy.testing as npt

TST_DATA = Path(__file__).parent / "data" / "test_data"


def test_read_spectral_bands():
    expected_spectral_band_index = [1, 2, 3, 4, 5]
    expected_lower_wavelength_limit = [
        0.0010000000474974513,
        0.0005000000237487257,
        0.00033333332976326346,
        0.0002500000118743628,
        0.00019999999494757503,
    ]
    expected_upper_wavelength_limit = [
        0.009999999776482582,
        0.0010000000474974513,
        0.0005000000237487257,
        0.00033333332976326346,
        0.0002500000118743628,
    ]
    actual = synthobs.read_spectral_bands(TST_DATA / "spectral" / "sp_sw_500ir_bd_hatp11")
    assert isinstance(actual, np.ndarray)
    assert actual.shape[0] == 500
    npt.assert_allclose(expected_spectral_band_index, actual["spectral_band_index"][:5])
    npt.assert_allclose(expected_lower_wavelength_limit, actual["lower_wavelength_limit"][:5])
    npt.assert_allclose(expected_upper_wavelength_limit, actual["upper_wavelength_limit"][:5])


def test_read_normalized_stellar_flux():
    expected = [
        3.7919454902446414e-12,
        7.004032465118826e-09,
        1.939694804775627e-08,
        3.824304428690084e-08,
        6.35967865036946e-08,
    ]
    actual = synthobs.read_normalized_stellar_flux(TST_DATA / "spectral" / "sp_sw_500ir_bd_hatp11")
    assert isinstance(actual, iris.cube.Cube)
    assert actual.shape[0] == 500
    npt.assert_allclose(expected, actual.data[:5])


def test_calc_stellar_flux():
    expected = [4.66409295e-10, 8.61495993e-07, 2.38582461e-06, 4.70389445e-06, 7.82240474e-06]
    actual = synthobs.calc_stellar_flux(TST_DATA / "spectral" / "sp_sw_500ir_bd_hatp11", 123)
    assert isinstance(actual, iris.cube.Cube)
    assert actual.shape[0] == 500
    npt.assert_allclose(expected, actual.data[:5])


def test_calc_transmission_spectrum_day_night_average():
    expected_spectral_band_centers = [
        0.0055,
        0.00075,
        0.00041666668,
        0.00029166666,
        0.000225,
    ]
    expected_rp_eff_over_rs = [0.99991664, 0.99993102, 0.99993353, 0.99994311, 0.99994491]
    (
        actual_spectral_band_centers,
        actual_rp_eff_over_rs,
    ) = synthobs.calc_transmission_spectrum_day_night_average(
        TST_DATA / "spectral" / "sp_sw_500ir_bd_hatp11",
        123,
        123,
        123,
        iris.load_cube(str(TST_DATA / "netcdf" / "planet_transmission_day.nc")),
        iris.load_cube(str(TST_DATA / "netcdf" / "planet_transmission_night.nc")),
    )
    assert isinstance(actual_spectral_band_centers, np.ndarray)
    assert isinstance(actual_rp_eff_over_rs, iris.cube.Cube)
    assert actual_spectral_band_centers.shape[0] == 500
    assert actual_rp_eff_over_rs.shape[0] == 500
    npt.assert_allclose(expected_spectral_band_centers, actual_spectral_band_centers[:5])
    npt.assert_allclose(expected_rp_eff_over_rs, actual_rp_eff_over_rs.data[:5])
