"""
setup.py:
    This file is used to install the package.
"""

import os
import glob
from setuptools import setup, find_packages, Extension

print("Running setup.py...")


# Set the project metadata
NAME = "contextsv"
VERSION = "1.0.0"
AUTHOR = "WGLab"
DESCRIPTION = "ContextSV: A tool for integrative structural variant detection."

# Get the conda environment's include path
conda_prefix = os.environ.get("CONDA_PREFIX")
if conda_prefix is None:
    raise AssertionError("CONDA_PREFIX is not set.")

conda_include_dir = os.path.join(conda_prefix, "include")

# Get the conda environment's lib path
conda_lib_dir = os.path.join(conda_prefix, "lib")

# Set the project dependencies
SRC_DIR = "src"
# SRC_FILES = glob.glob(os.path.join(SRC_DIR, "*.cpp"))
SRC_FILES = glob.glob(os.path.join(SRC_DIR, "*.cpp"))
SRC_FILES = [f for f in SRC_FILES if "main.cpp" not in f]  # Ignore the main.cpp file

INCLUDE_DIR = "include"
INCLUDE_FILES = glob.glob(os.path.join(INCLUDE_DIR, "*.h"))

# Set up the extension
ext = Extension(
    name="_" + NAME,
    sources=SRC_FILES,
    include_dirs=[INCLUDE_DIR, conda_include_dir],
    extra_compile_args=["-std=c++17"],
    language="c++",
    libraries=["hts"],
    library_dirs=[conda_lib_dir]
)

# Set up the module
setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    description=DESCRIPTION,
    ext_modules=[ext],
    py_modules=[NAME],
    packages=find_packages(),
    test_suite="tests",
    entry_points={
        "console_scripts": [
            "contextsv = contextsv:main"
        ]
    }
)

