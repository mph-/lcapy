"""This module defines special symbols f, s, t, and omega.

Copyright 2014--2021 Michael Hayes, UCECE

"""

# The following domain variables are Lcapy expressions.  For example, f is an
# FourierDomainExpression object wrapping the fsym symbol.
from .fexpr import f
from .texpr import t
from .sexpr import s
from .omegaexpr import omega
from .normomegaexpr import Omega
from .jfexpr import jf, j2pif
from .jomegaexpr import jomega
from .normfexpr import F
from .nexpr import n
from .kexpr import k
from .zexpr import z

# Import common SymPy symbols.
from .sym import pi, j, oo, inf, one, omega0sym, f0sym

# Perhaps have class for domainconstants such as omega0, t0, f0, etc.
omega0 = omega0sym
f0 = f0sym

jw = jomega

jomega0 = j * omega0
jw0 = jomega0

# This represents an arbitrary angular frequency
w = omega
# This represents a specific angular frequency and is assumed to be positive
w0 = omega0

W = Omega

domain_vars = [t, f, s, omega, omega0, Omega, F, n, k, z, jw, jw0]
domain_var_ids = [id(var) for var in domain_vars]
