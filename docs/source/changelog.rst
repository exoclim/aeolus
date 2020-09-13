Changelog
=========

.. default-role:: py:obj

0.4.6
-----

:Date: TBA

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
  - Refactor surface and TOA energy balance calculation, and do notapply spatial averaging to P-E
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
* API changes: merge `grid` into `coord`; move `misc` flux calculations to `calc` submodule,
remove `util` folder by moving `text` to `plot`.
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
