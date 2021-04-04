#CYTHON -- MAC OS X FIX --	following https://github.com/cython/cython/issues/1725
import numpy
import os
import pyximport
numpy_path = numpy.get_include()
os.environ['CFLAGS'] = "-I" + numpy_path
pyximport.install(
	language_level=3,
    pyimport=False,
    setup_args={'include_dirs': numpy.get_include()}
    )

from . import Cfourvec as Cfv

# Definition modules
from . import const
from . import fourvec

