#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Package build and install script."""
from setuptools import setup

import versioneer


CMDCLASS = versioneer.get_cmdclass()


setup(
    version=versioneer.get_version(),
    cmdclass=CMDCLASS,
)
