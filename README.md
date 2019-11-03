[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![Conda](https://img.shields.io/conda/v/dennissergeev/aeolus?color=dark-green&logo=anaconda)](https://anaconda.org/dennissergeev/aeolus)
[![PyPI](https://img.shields.io/pypi/v/aeolus.svg?logo=python&logoColor=white)](https://pypi.org/project/aeolus/)
[![Build Status](https://travis-ci.com/exoclim/aeolus.svg?branch=master)](https://travis-ci.com/exoclim/aeolus?branch=master)
[![Documentation Status](https://readthedocs.org/projects/aeolus/badge/?version=latest)](https://aeolus.readthedocs.io/en/latest/?badge=latest)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](LICENSE)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# aeolus :wind_face:
Python library for object-oriented analysis of atmospheric model output built on top of [iris](https://github.com/SciTools/iris).

The documentation is available [here](https://aeolus.readthedocs.io/en/latest/).

Contributions are welcome.

## Installation

The list of dependencies can be found in [`ci/`](ci/requirements-py37.yml) and the best way to install them in a separate environment is to use [`conda`](https://conda.io)
```bash
conda install -c conda-forge cartopy iris matplotlib numpy
pip install latlon23
```

After the required packages are installed, `aeolus` can be installed either from PyPI or from source.

### PyPI
```bash
pip install aeolus
```


### From source
```bash
git clone https://github.com/exoclim/aeolus.git

cd aeolus

python setup.py install
```
