"""This module provides current support.

Copyright 2020 Michael Hayes, UCECE

"""
from .expr import expr
from .sym import omega0sym
from .symbols import s, omega0
from .cexpr import ConstantExpression, ConstantCurrent
from .fexpr import FourierDomainExpression, FourierDomainCurrent
from .omegaexpr import AngularFourierDomainExpression, AngularFourierDomainCurrent
from .sexpr import LaplaceDomainExpression, LaplaceDomainCurrent
from .texpr import TimeDomainExpression, TimeDomainCurrent
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


def current(arg, **assumptions):

    mapping = {ConstantExpression: ConstantCurrent,
               TimeDomainExpression: TimeDomainCurrent,
               LaplaceDomainExpression: LaplaceDomainCurrent,
               FourierDomainExpression: FourierDomainCurrent,
               AngularFourierDomainExpression: AngularFourierDomainCurrent}
    
    expr1 = expr(arg, **assumptions)
    if expr1.__class__ in mapping:
        expr1 = mapping[expr1.__class__](expr1)

    return expr1.apply_unit(uu.amperes)
