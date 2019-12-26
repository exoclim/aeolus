#!/usr/bin/env python
"""Package build and install script."""
import os
import sys

from setuptools import find_packages, setup

import versioneer


# Publish the library to PyPI.
if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist upload")
    sys.exit()


def get_readme():
    """Load README for display on PyPI."""
    with open("README.md") as f:
        return f.read()


CMDCLASS = versioneer.get_cmdclass()


setup(
    name="aeolus",
    version=versioneer.get_version(),
    cmdclass=CMDCLASS,
    description="Python library for object-oriented analysis of atmospheric model output",
    long_description=get_readme(),
    long_description_content_type="text/markdown",
    author="Denis Sergeev",
    author_email="dennis.sergeev@gmail.com",
    url="https://github.com/exoclim/aeolus",
    package_dir={"aeolus": "aeolus"},
    include_package_data=True,
    # package_data={"aeolus": ["aeolus/const/store/*.json", "aeolus/tests/data/*.json"]},
    packages=find_packages(),
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.16",
        "pytest>=3.3",
        "matplotlib>=2",
        "scitools-iris>=2.3",
        "latlon23",
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
)
