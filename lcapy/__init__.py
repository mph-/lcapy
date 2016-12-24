"""
Lcapy is a Python library for linear circuit analysis.  It uses SymPy
for symbolic mathematics.

Lcapy can analyse circuits described with netlists using modified nodal
analysis.  See lcapy.netlist

Alternatively, Lcapy can analyse networks and circuits formed by
combining one, two, and three port networks.

Copyright 2014, 2015, 2016 Michael Hayes, UCECE
"""

from __future__ import absolute_import, print_function
from pkg_resources import get_distribution

__version__ = get_distribution('lcapy').version

import sys
if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    raise ImportError("Python Version 2.6 or above is required for Lcapy.")
else:  # Python 3
    pass
    # Here we can also check for specific Python 3 versions, if needed

del sys

from sympy import init_printing
init_printing()


# List of symbols that get imported with 'from lcapy import *'
__all__ = []

# Add the modules that get searched to allow 'from lcapy import V'
# rather then having to specify module, 'from lcapy.oneport import V'
from . import core
__all__.extend(core.__all__)
from .core import *

from . import oneport
__all__.extend(oneport.__all__)
from .oneport import *

#from . import twoport
#__all__.extend(twoport.__all__)
#from .twoport import *

#from . import threeport
#__all__.extend(threeport.__all__)
#from .threeport import *

from . import circuit
__all__.extend(circuit.__all__)
from .circuit import *

from . import schematic
__all__.extend(schematic.__all__)
from .schematic import *

__all__.extend(('show_version', ))

def show_version():
    
    from sys import version as python_version
    from sympy import __version__ as sympy_version
    from numpy import __version__ as numpy_version

    print('Python: %s\nSymPy: %s\nNumPy: %s' % 
          (python_version, sympy_version, numpy_version))

