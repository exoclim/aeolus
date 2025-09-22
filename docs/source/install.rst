.. _install:

############
Installation
############

.. note::

The main dependency of aeolus is `iris <https://scitools.org.uk/iris/docs/latest/>`_, but various
minor parts also depend on the following Python packages

- cached-property
- cartopy
- matplotlib
- numpy
- latlon23
- metpy
- python-stratify
- pyvista
- xarray


Recommended method: pixi
========================

Install `pixi <https://pixi.sh/>`_ and create a local project following the instructions there.
Add aeolus to the project:

  pixi add aeolus

Verify the installation:

  pixi run aeolus --version


Alternative method 1: conda
===========================

Another installation method is via `conda <https://conda.io/docs/>`_ ::

  conda install -c conda-forge aeolus

Alternative method 2: pip
=========================
You can also install aeolus from the Python Package Index (PyPI) ::

  pip install aeolus


Alternative method 3: install from source
=========================================

To get the latest (potentially unstable) version of the library you can directly clone the `GitHub repository <https://github.com/exoclim/aeolus>`_ ::

  git clone https://www.github.com/exoclim/aeolus.git
  cd aeolus

and install aeolus in the standard mode ::

  python setup.py install

and install aeolus in the developer mode ::

  python setup.py develop

or::

  pip install -e .


Verifying installation
======================

Once installed via any of these methods, you can run aeolus's suite of
tests using `pytest <http://doc.pytest.org/>`_.  From the top-level
directory of the aeolus installation ::

  conda install pytest  # if you don't have it already; or 'pip install pytest', or 'pixi add pytest'
  pytest aeolus

If you don't know the directory where aeolus was installed, you can find it via ::

  python -c "import aeolus; print(aeolus.__path__[0])"

If the pytest command results in any error messages or test failures,
something has gone wrong, and please refer to the Troubleshooting
information below.

Troubleshooting
===============

Please search through the `Issues page`_ on Github if anybody else has had the same problem you're facing.
If none do, then please send open a new Issue.

.. _Issues page: https://github.com/exoclim/aeolus/issues
