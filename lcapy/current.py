"""This module provides current support.

Copyright 2020 Michael Hayes, UCECE

"""
from .expr import expr
from .sym import omega0sym
from .symbols import s, omega0
from .cexpr import ConstantCurrent
from .fexpr import FourierDomainCurrent
from .omegaexpr import AngularFourierDomainCurrent
from .sexpr import LaplaceDomainCurrent
from .texpr import TimeDomainCurrent
from .noiseomegaexpr import AngularFourierDomainNoiseCurrent
from .noisefexpr import FourierDomainNoiseCurrent
from .phasor import PhasorCurrent
from .units import u as uu


def Iname(name, kind, cache=False):

    if kind in ('s', 'laplace'):    
        return LaplaceDomainCurrent(name + '(s)')
    elif kind in ('t', 'time'):
        return TimeDomainCurrent(name.lower() + '(t)')
    elif kind in (omega0sym, omega0, 'ac'):
        return PhasorCurrent(name + '(omega_0)')
    # Not caching is a hack to avoid conflicts of In1 with In1(s) etc.
    # when using subnetlists.  The alternative is a proper context
    # switch.  This would require every method to set the context.
    return expr(name, cache=cache)            


def Itype(kind):
    
    if isinstance(kind, str) and kind[0] == 'n':
        return AngularFourierDomainNoiseCurrent
    try:
        return {'ivp' : LaplaceDomainCurrent,
                's' : LaplaceDomainCurrent,
                'n' : AngularFourierDomainNoiseCurrent,
                'ac' : PhasorCurrent,
                'dc' : ConstantCurrent,
                't' : TimeDomainCurrent,
                'time' : TimeDomainCurrent}[kind]
    except KeyError:
        return PhasorCurrent


def current(arg):

    expr1 = expr(arg)
    return expr1.apply_unit(uu.amperes)
