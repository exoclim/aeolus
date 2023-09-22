Changelog
=========

.. default-role:: py:obj


0.4.18
------

:Date: 22 September 2023

* Re-fix the bug in day-night transmission spectra averaging.
* Add `lfric` functions in a submodule with the same name.
* Add `esmf_regrid` as a dependency.
* Move `calc_derived_cubes` from `pouch`, add `air_pressure()`.
* Add various plotting utilities.
* Add the `log` submodule with `loguru` as a dependency.
* Add `proc_um_output`.
* Add the `RUNTIME` variable.
* Update dependencies.


0.4.17
------

:Date: 14 April 2023

* Add `lfric` name dictionary.
* Add `subset.CellMethodConstraint`.
* Update `pyvista` functions to be in step with the latest `pyvista` version.
* Fix a bug in the `calc_transmission_spectrum()` (does not affect the results).


0.4.16
------

:Date: 11 January 2023

* Miscellaneous improvements.


0.4.15
------

:Date: 22 April 2022

* Patch the calculation of the day-night average transmission flux.
* Update citation.


0.4.14
------

:Date: 30 March 2022

* Require iris 3.2 as a minimum.
* Drop Python 3.7 support.
* Add github action tests for Python 3.10 (for now excluding windows).
* Add new model names.
* Correct the calculation of `sigma_p` and remove its averaging.
* Make `calc.stats.time_mean()` be able to accept cubes or collections of cubes.
* Add the `structured_um_loading` option to `io.load_data()`.
* Add `wind_rot_div()` to the diagnostics module and a relevant example notebook.
* Simplify `isel()` following suggestions from the `iris` devs.


0.4.13
------

:Date: 18 October 2021

* Make the package compatible with iris v3.1 by using explicit imports
* Refactor functions to calculate transmission spectrum.
* Add `io.create_dummy_cube()`
* Add `synthobs.calc_transmission_spectrum()`
* Use `da.roll()` instead of `np.roll` in `coord.roll_*` functions.
* Fix  a few minor bugs in synthobs.


0.4.12
------

:Date: 29 July 2021

* Add `coord_check` option to `AtmoSim` and `CoordContainer` to allow for non-checking of coordinates.
* Remove `MidpointNormalize`
* Add `plot.text.tex2cf_units()` and `plot.text.unit_format()`
* Rewrite github actions
* Add pre-commit hooks


0.4.11
------

:Date: 08 June 2021

* Add representations of `AtmoSim`
* Fix a bug with the `const_from_attrs` decorator
* Add contributor's guide
* Replace all autodoc function and class references by module references


0.4.10
------

:Date: 07 June 2021

* Add functions to calculate synthetic transmission spectrum. By `Maria Zamyatina <https://github.com/mzamyatina>`_.
* Replace `interp_to_pres_lev()`, `interp_all_to_pres_lev()` and `interp_to_single_pres_lev()` with `interp_cube_from_height_to_pressure_levels()` and `interp_cubelist_from_height_to_pressure_levels()` with a better interface
* Move test data to a separate repository: `aeolus_data`
* Add more names to `model.base`
* Add an option for the `const_from_attrs()` decorator not to raise an error


0.4.9
-----

:Date: 12 April 2021

* Add cached-property as a dependency
* Rename `AtmosFlow` to `AtmoSim` and create a base class `AtmoSimBase`
* Add `extract()` method to `AtmoSimBase`
* Add pressure coordinate to `DimConstr`
* Refactor `Run` and prepare for its deprecation
* Add `load_data()` to `io`
* Move `add_planet_conf_to_cubes()` to the `const` module
* Deprecate `ScalarCube`
* Add new variable names to `um`
* Refactor derived constants and add `planet_rotation_rate` to the recipes
* Add an option not to broadcast the coordinate to the cube's shape in `coord_to_cube()`
* Make `spatial()`, `time_mean()` and `vertical_mean()` return the input cube in case of `CoordinateCollapseError`
* Add `abs_coord_mean()` to average data over latitudes symmetric around the equator
* Add functions to calculate meridional and zonal streamfunctions
* Improve docstrings
* Add an example notebook for working with model names


0.4.8
-----

:Date: 25 January 2021

* Adapt to iris v3.0
* Add new meta decorators
* Fix typos


0.4.7
-----

:Date: 03 December 2020

* Move to conda-forge for building the package

* Replace TravisCI with GitHub Actions

* Restructure the package:

  - the library is now in `src/aeolus`
  - tests are now in `tests/`

* Core classes:

  - Add `AtmosFlow`

* New `calc` functions:

  - Add a decorator to update cube metadata, `update_metadata()`
  - Add shortcut functions `spatial_mean()` and `time_mean()`
  - Add `air_potential_temperature()`
  - Add functions to rotate and regrid variables to "tidally-locked" coordinates
  - Add `wind_speed()`

* New `coord` functions:

  - Add functions to interpolate cubes to pressure levels (depend on python-stratify package)
  - Add a function to interpolate one cube to another along the time dimension (`interp_to_cube_time`)
  - Add a container to store common coordinates
  - Add `check_coords()`, `get_xy_coords()`

* New `subset` functions:

  - Add a function to filter out duplicated cubes from a cubelist: `unique_cubes()`

* Other changes:

  - Append names to `model.um`
  - Add a function to load vertical levels data
  - Improve `interp_to_pres_lev()`
  - Add `model` keyword to `plot.pv` functions
  - Rewrite `DimConstr` API


0.4.6
-----

:Date: 17 September 2020

* New `calc` functions:

  - Add convenience functions to `calculus` for doing simple derivatives
  - Add a function to calculate horizontal divergence in spherical coordinates
  - Add `air_density()`, `air_temperature()`, `dry_lapse_rate()`, `flux()` and `geopotential_height` functions
  - Add `horiz_wind_cmpnts()` helper function
  - Add `normalize_cube()`
  - Add `superrotation_index()`

* New `coord` functions:

  - Add a function to emulate `xarray`'s `isel()` method.
  - Update `get_cube_datetimes()` and add a new function, `get_cube_rel_days()`
  - Add a function to broadcast coordinate deltas to a cube.
  - Add a function to calculate volume from a cube's grid.

* Other changes:

  - Append names to `model.um`
  - Refactor surface and TOA energy balance calculation, and do not apply spatial averaging to P-E
  - Override `__repr__` of `model.base.Model`

0.4.5
-----

:Date: 08 June 2020

* API changes:

  - add `model` submodule for model-specific variable and coordinate names
  - replace all `UM_*` variables with the `model` reference
  - replace all `DIM_CONSTR_*` by a class `DimConstr` with each of the constraints as an attribute
  - replace `coord.add_binned_lon_lat()` by a generic `coord.add_binned_coord()` function

* Minor bug fixes and clean-up

0.4.4
-----

:Date: 04 May 2020

* Add Python 3.8 to build matrix
* API changes: merge `grid` into `coord`; move `misc` flux calculations to `calc` submodule, remove `util` folder by moving `text` to `plot`.
* Add a function to attach non-Earth auxiliary time coordinates to a cube
* Add Titan constants (some orbital parameters are those for Saturn for simplicity)
* Fix a typo in Earth constants
* Add a new function for matplotlib plots: `plot.add_custom_legend()`
* Minor fixes in the travis integration

0.4.3
-----

:Date: 30 March 2020

* Add diagnostics: `vertical_mean()`, `vertical_sum()`
* Improve diagnostics (`sfc_water_balance()`) and utilities (`regrid3d()`)
* Allow for the initialisation of `Run` from a pre-processed data
* Add a method to `Run` to save processed cubelist to netCDF
* Fix a few bugs

0.4.2
-----

:Date: 05 January 2020

* Improve calculation of precipitation sums
* Add a helper function to retrieve planet radius from a cube
* Add `timestep` attribute to `Run`
* Improve docstrings
* Remove two functions from `pv` submodule (now in `pyvista` library)
* Move documentation to github pages
* Fix a few bugs

0.4.1
-----

:Date: 03 December 2019

* Add a few standard constants
* Improve units in metpy interface


0.4.0
-----

:Date: 28 November 2019

* Add metpy-to-iris interface
* Fix a few bugs


0.3.2
-----

:Date: 21 November 2019

* Add basic examples as Jupyter Notebooks
* Improve plotting functions and diags
* Add test data


0.2
---

:Date: 02 November 2019

* Technical updates

0.1
---

:Date: 31 October 2019

* First packaged release
