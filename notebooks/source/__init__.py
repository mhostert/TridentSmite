#CYTHON -- MAC OS X FIX -- https://github.com/cython/cython/issues/1725
import numpy as np
import os
import pyximport

numpy_path = np.get_include()
os.environ['CFLAGS'] = "-I" + numpy_path
pyximport.install(
	language_level=3,
    pyimport=False,
    setup_args={'include_dirs': np.get_include()}
    )

from . import Cfourvec as Cfv
from . import const
from . import kinematics as kin
from . import parser
from . import generator as gen
from . import exp_classes as exps
