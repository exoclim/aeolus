# -*- coding: utf-8 -*-
"""Test synthobs submodule."""
from pathlib import Path

from aeolus import synthobs

import iris
import iris.coords
import iris.cube

import numpy as np
import numpy.testing as npt

import pytest


iris.FUTURE.datum_support = True
TST_DATA = Path(__file__).parent / "data" / "test_data"


@pytest.fixture(scope="module")
def example_trans_day():
    return iris.load_cube(str(TST_DATA / "netcdf" / "planet_transmission_day.nc"))


@pytest.fixture(scope="module")
def example_trans_night():
    return iris.load_cube(str(TST_DATA / "netcdf" / "planet_transmission_night.nc"))


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


def test_calc_geom_mean_mirrored(example_trans_day, example_trans_night):
    actual = synthobs.calc_geom_mean_mirrored(example_trans_day, example_trans_night)
    npt.assert_allclose(actual.data.max(), 1.7307187876145394e-08)
    npt.assert_allclose(
        actual.data[123, 1, 71:74],
        [2.40508136e-11, 2.40748625e-11, 2.40532320e-11],
    )
    assert actual.units == example_trans_day.units
    assert actual.shape == example_trans_day.shape


def test_calc_transmission_spectrum(example_trans_day):
    expected_spectral_bands = [
        0.0055,
        0.00075,
        0.00041666668,
        0.00029166666,
        0.000225,
    ]
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
    expected_rp_eff_over_rs = [0.06043149, 0.0605353, 0.06054769, 0.060621, 0.06062967]
    actual_rp_eff_over_rs = synthobs.calc_transmission_spectrum(
        example_trans_day,
        TST_DATA / "spectral" / "sp_sw_500ir_bd_hatp11",
        300.85538168664425,
        475026500.0,
        28995379.0,
    )
    actual_spectral_bands = actual_rp_eff_over_rs.coord("spectral_band_centres")
    actual_spectral_band_ll = actual_rp_eff_over_rs.coord("spectral_band_lower_limit")
    actual_spectral_band_ul = actual_rp_eff_over_rs.coord("spectral_band_upper_limit")
    assert isinstance(actual_spectral_bands, iris.coords.AuxCoord)
    assert isinstance(actual_rp_eff_over_rs, iris.cube.Cube)
    assert actual_spectral_bands.shape[0] == 500
    assert actual_rp_eff_over_rs.shape[0] == 500
    npt.assert_allclose(expected_spectral_bands, actual_spectral_bands.points[:5])
    npt.assert_allclose(expected_lower_wavelength_limit, actual_spectral_band_ll.points[:5])
    npt.assert_allclose(expected_upper_wavelength_limit, actual_spectral_band_ul.points[:5])
    npt.assert_allclose(expected_rp_eff_over_rs, actual_rp_eff_over_rs.data[:5])


def test_calc_transmission_spectrum_day_night_average(example_trans_day, example_trans_night):
    expected_rp_eff_over_rs = [
        0.9999097552081786,
        0.9999249847710132,
        0.9999267911653156,
        0.9999376935302737,
        0.9999389664305547,
    ]
    actual_rp_eff_over_rs = synthobs.calc_transmission_spectrum_day_night_average(
        example_trans_day,
        example_trans_night,
        TST_DATA / "spectral" / "sp_sw_500ir_bd_hatp11",
        123,
        123,
        123,
    )
    assert isinstance(actual_rp_eff_over_rs, iris.cube.Cube)
    assert actual_rp_eff_over_rs.shape[0] == 500
    npt.assert_allclose(expected_rp_eff_over_rs, actual_rp_eff_over_rs.data[:5])
