"""
Lcapy is a Python library for symbolic linear circuit and signal analysis.

Lcapy can analyse circuits described with netlists using modified nodal
analysis.  See lcapy.netlist

Alternatively, Lcapy can analyse networks and circuits formed by
combining one, two, and three port networks.  See lcapy.oneport

For detailed documentation see http://lcapy.readthedocs.io/en/latest

Copyright 2014--2022 Michael Hayes, UCECE
"""
from __future__ import absolute_import, print_function
from sympy import Symbol
from sympy.core.sympify import converter
# This must be imported early to avoid circular import with expr
from .functions import *
from .units import volts, amperes, ohms, siemens, watts
from .state import state
from .inverse_dft import *
from .dft import *
from .differenceequation import *
from .dltifilter import *
from .differentialequation import *
from .ltifilter import *
from .discretetime import *
from .phasor import phasor, phasor_ratio
from .normfexpr import Fexpr
from .normomegaexpr import Omegaexpr
from .omegaexpr import omegaexpr
from .cexpr import cexpr
from .texpr import texpr
from .sexpr import sexpr, zp2tf, tf, pr2tf
from .fexpr import fexpr
from .expr import *
from .simulator import *
from .randomnetwork import *
from .nettransform import *
from .laplace import *
from .dtstatespace import *
from .statespace import *
from .vector import *
from .tmatrix import *
from .smatrix import *
from .matrix import *
from .sym import *
from .printing import *
from .susceptance import susceptance
from .reactance import reactance
from .inductance import inductance
from .capacitance import capacitance
from .conductance import conductance
from .resistance import resistance
from .transfer import transfer
from .impedance import impedance
from .admittance import admittance
from .current import current, noisecurrent, phasorcurrent
from .voltage import voltage, noisevoltage, phasorvoltage
from .twoport import *
from .oneport import *
from .circuit import *
from .symbols import *
import sys
import pkg_resources
del absolute_import, print_function

name = "lcapy"


__version__ = pkg_resources.require('lcapy')[0].version
lcapy_version = __version__

if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    raise ImportError("Python Version 2.6 or above is required for Lcapy.")
else:  # Python 3
    pass
    # Here we can also check for specific Python 3 versions, if needed

del sys

# Do not import units.u since this will conflict with unit step


def show_version():
    """Show versions of Lcapy, SymPy, NumPy, MatplotLib, SciPy, and Python."""

    from sys import version as python_version
    from sympy import __version__ as sympy_version
    from numpy import __version__ as numpy_version
    from scipy import __version__ as scipy_version
    from matplotlib import __version__ as matplotlib_version

    print('Python: %s\nSymPy: %s\nNumPy: %s\nMatplotlib: %s\nSciPy: %s\nLcapy: %s' %
          (python_version, sympy_version, numpy_version,
           matplotlib_version, scipy_version, lcapy_version))


# The following is to help sympify deal with j.
# A better fix might be to define an Lcapy class for j and to
# use the __sympy_ method.
converter['j'] = j
converter[Symbol('j')] = j
del converter, Symbol
