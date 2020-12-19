"""This module provides voltage support.

Copyright 2020 Michael Hayes, UCECE

"""
from .expr import expr
from .sym import omega0sym
from .symbols import s, omega0
from .cexpr import ConstantExpression, ConstantVoltage
from .fexpr import FourierDomainExpression, FourierDomainVoltage
from .omegaexpr import AngularFourierDomainExpression, AngularFourierDomainVoltage
from .sexpr import LaplaceDomainExpression, LaplaceDomainVoltage
from .texpr import TimeDomainExpression, TimeDomainVoltage
from .noiseomegaexpr import AngularFourierDomainNoiseVoltage
from .noisefexpr import FourierDomainNoiseVoltage
from .phasor import PhasorVoltage
from .units import u as uu


def Vname(name, kind, cache=False):

    if kind in ('s', 'laplace'):    
        return LaplaceDomainVoltage(name + '(s)')
    elif kind in ('t', 'time'):
        return TimeDomainVoltage(name.lower() + '(t)')
    elif kind in (omega0sym, omega0, 'ac'):
        return PhasorVoltage(name + '(omega_0)')
    # Not caching is a hack to avoid conflicts of Vn1 with Vn1(s) etc.
    # when using subnetlists.  The alternative is a proper context
    # switch.  This would require every method to set the context.
    return expr(name, cache=cache)            


def Vtype(kind):
    
    if isinstance(kind, str) and kind[0] == 'n':
        return AngularFourierDomainNoiseVoltage
    try:
        return {'ivp' : LaplaceDomainVoltage,
                's' : LaplaceDomainVoltage,
                'n' : AngularFourierDomainNoiseVoltage,
                'ac' : PhasorVoltage,
                'dc' : ConstantVoltage,
                't' : TimeDomainVoltage,
                'time' : TimeDomainVoltage}[kind]
    except KeyError:
        return PhasorVoltage


def voltage(arg):

    mapping = {ConstantExpression: ConstantVoltage,
               TimeDomainExpression: TimeDomainVoltage,
               LaplaceDomainExpression: LaplaceDomainVoltage,
               FourierDomainExpression: FourierDomainVoltage,
               AngularFourierDomainExpression: AngularFourierDomainVoltage}
    
    expr1 = expr(arg)
    if expr1.__class__ in mapping:
        expr1 = mapping(expr1)
    
    return expr1.apply_unit(uu.volts)
