#! /usr/bin/env python3

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy
import os

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

# extension module
_fast_mutation = Extension("_fast_mutation",
                   ["wrapper.i","fast_mutation.cpp"],
                   language = "c++",
                   extra_compile_args=["-std=c++11", "-O3"],
                   include_dirs = [numpy_include],
                   )

# setup
setup(  name        = "fast_mutation module",
        description = "a fast version of mutation and expand",
        author      = "Bohao Tang",
        version     = "1.0",
        ext_modules = [_fast_mutation]
        )
