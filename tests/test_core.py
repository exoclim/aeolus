# -*- coding: utf-8 -*-
"""Test the core submodule."""
# import contextlib
# import itertools
# import shutil
# import tempfile
# from datetime import datetime
from pathlib import Path

#
# from aeolus import core
from aeolus.exceptions import LoadError

# import iris
#
# import numpy as np
# import numpy.testing as npt

import pytest

TST_DATA = Path(__file__).parent / "data"

# _counter = itertools.count()


def test_foo():
    pass


def test_loaderror():
    with pytest.raises(LoadError):
        raise LoadError
