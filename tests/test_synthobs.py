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
    return iris.load_cube(TST_DATA / "netcdf" / "planet_transmission_day.nc")


@pytest.fixture(scope="module")
def example_trans_night():
    return iris.load_cube(TST_DATA / "netcdf" / "planet_transmission_night.nc")


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
    # Default
    actual = synthobs.calc_geom_mean_mirrored(example_trans_day, example_trans_night)
    npt.assert_allclose(actual.data.max(), 1.5901232591606386e-08)
    npt.assert_allclose(
        actual.data[123, 45, 102:108],
        [
            0.00000000e00,
            5.28650299e-40,
            1.09253009e-10,
            7.16136321e-10,
            5.28943014e-10,
            0.00000000e00,
        ],
    )
    assert actual.units == example_trans_day.units
    assert actual.shape == example_trans_day.shape
    # With an additional shift along the x-coordinate
    actual = synthobs.calc_geom_mean_mirrored(example_trans_day, example_trans_night, add_shift=-1)
    npt.assert_allclose(actual.data.max(), 1.76899099563551e-08)
    npt.assert_allclose(
        actual.data[123, 45, 102:109],
        [
            0.00000000e00,
            1.54026439e-67,
            2.12892023e-11,
            6.61482043e-10,
            7.37627923e-10,
            3.77198057e-10,
            0.00000000e00,
        ],
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
    # Default
    expected_rp_eff_over_rs = [
        0.1289366838346542,
        0.1289351109234838,
        0.1289354291294527,
        0.1289365978678745,
        0.1289363841112755,
    ]
    actual_rp_eff_over_rs = synthobs.calc_transmission_spectrum_day_night_average(
        example_trans_day,
        example_trans_night,
        spectral_file=TST_DATA / "spectral" / "sp_sw_500ir_bd_hatp11",
        stellar_constant_at_1_au=1272.86475403,
        stellar_radius=7.302834e08,
        planet_top_of_atmosphere=94193200,
    )
    assert isinstance(actual_rp_eff_over_rs, iris.cube.Cube)
    assert actual_rp_eff_over_rs.shape[0] == 500
    npt.assert_allclose(expected_rp_eff_over_rs, actual_rp_eff_over_rs.data[-5:])
    # With an additional shift
    expected_rp_eff_over_rs = [
        0.1289324245796913,
        0.1289306558725381,
        0.1289309758113692,
        0.1289323585757668,
        0.1289320864595353,
    ]
    actual_rp_eff_over_rs = synthobs.calc_transmission_spectrum_day_night_average(
        example_trans_day,
        example_trans_night,
        add_shift=-1,
        spectral_file=TST_DATA / "spectral" / "sp_sw_500ir_bd_hatp11",
        stellar_constant_at_1_au=1272.86475403,
        stellar_radius=7.302834e08,
        planet_top_of_atmosphere=94193200,
    )
    npt.assert_allclose(expected_rp_eff_over_rs, actual_rp_eff_over_rs.data[-5:])
