.. _install:

############
Installation
############

Supported language:

- Python: 3.7

.. note::

   We highly recommend installing and using the free `Anaconda
   <https://www.anaconda.com/download/>`_ distribution of Python (or
   `Miniconda <https://conda.io/miniconda.html>`_, if you don't want
   all of the extra packages that come built-in with Anaconda), which
   works on Mac, Linux, and Windows, both on normal computers and
   institutional clusters and doesn't require root permissions.

Aeolus depends on the following Python packages

- cartopy
- iris
- matplotlib
- numpy
- latlon23
- metpy
- xarray

Additionally, a :ref:`submodule <pyvista_ref>` that provides an interface to `PyVista <https://docs.pyvista.org/>`_ obviously requires it to be installed.

After the required packages are installed, aeolus can be installed either from Anaconda, PyPI, or from source.


Recommended installation method: conda
======================================

The recommended installation method is via `conda <https://conda.io/docs/>`_

First, prepare the environment by conda-installing all the dependencies from the conda-forge channel ::

  conda install -c conda-forge cartopy iris matplotlib numpy metpy xarray
  pip install latlon23

Then install aeolus ::

  conda install -c dennissergeev aeolus

Alternative method: PyPI
========================
Install aeolus from the Python Package Index ::

  pip install aeolus


Alternative method: clone from Github
=====================================

You can also directly clone the `Github repo <https://github.com/exoclim/aeolus>`_ ::

  git clone https://www.github.com/exoclim/aeolus.git
  cd aeolus
  python setup.py install

Verifying proper installation
=============================

Once installed via any of these methods, you can run aeolus's suite of
tests using `py.test <http://doc.pytest.org/>`_.  From the top-level
directory of the aeolus installation ::

  conda install pytest  # if you don't have it already; or 'pip install pytest'
  py.test aeolus

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
