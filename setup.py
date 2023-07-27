"""
setup.py:
    This file is used to install the package.
"""

print("Running setup.py...")

import os
import glob
from setuptools import setup, find_packages, Extension


# Set the project metadata
name = "contextsv"
version = "0.0.1"
author = "WGLab"
description = "ContextSV: A tool for integrative structural variant detection."

# Get the conda environment's include path
conda_prefix = os.environ.get("CONDA_PREFIX")
if conda_prefix is None:
    raise Exception("CONDA_PREFIX is not set.")
conda_include_dir = os.path.join(conda_prefix, "include")

print("CONDA_PREFIX: {}".format(conda_prefix))  # DEBUG
print("include_dir: {}".format(conda_include_dir))  # DEBUG

# Get the conda environment's lib path
conda_lib_dir = os.path.join(conda_prefix, "lib")

print("lib_dir: {}".format(conda_lib_dir))  # DEBUG

# Set the project dependencies
src_dir = "src"
src_files = glob.glob(os.path.join(src_dir, "*.cpp"))
include_dir = "include"
include_files = glob.glob(os.path.join(include_dir, "*.h"))

# Set up the extension
ext = Extension(
    name="_" + name,
    sources=src_files,
    include_dirs=[include_dir, conda_include_dir],
    extra_compile_args=["-std=c++11"],
    language="c++",
    libraries=["hts"],
    library_dirs=[conda_lib_dir]
)

# Set up the module
setup(
    name=name,
    version=version,
    author=author,
    description=description,
    ext_modules=[ext],
    py_modules=[name],
    packages=find_packages(),
    test_suite="tests",
    entry_points={
        "console_scripts": [
            "contextsv = contextsv:main"
        ]
    }
)

