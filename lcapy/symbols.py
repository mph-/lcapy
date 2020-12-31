"""This module defines special symbols f, s, t, and omega.

Copyright 2014--2020 Michael Hayes, UCECE

"""

# The following are Lcapy expressions.  For example, f is an
# FourierDomainExpression object wrapping the fsym symbol.
from .fexpr import f
from .texpr import t
from .sexpr import s
from .omegaexpr import omega, omega0

# Import common SymPy symbols.
from .sym import pi, j, oo, inf, one

jomega = j * omega
jw = jomega

jomega0 = j * omega0
jw0 = jomega0
