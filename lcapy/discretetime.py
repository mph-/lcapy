"""This module provides discrete-time support.

It introduces three special variables:
   n for discrete-time sequences
   k for discrete-frequency sequences
   z for z-transforms.

Copyright 2020--2021 Michael Hayes, UCECE

"""

import sympy as sym
from .sym import sympify
from .nexpr import nexpr, n
from .kexpr import kexpr, k
from .zexpr import zexpr, z
from .dsym import nsym, ksym, zsym, dt, df

from .transform import transform as transform1
from .transform import call as call1
from .functions import Function
from .ztransform import *
from .seq import seq
