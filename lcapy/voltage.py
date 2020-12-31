"""This module provides voltage support.

Copyright 2020 Michael Hayes, UCECE

"""
from .expr import expr
from .sym import omega0sym
from .symbols import s, omega0
from .cexpr import ConstantVoltage
from .fexpr import FourierDomainVoltage
from .omegaexpr import AngularFourierDomainVoltage
from .sexpr import LaplaceDomainVoltage
from .texpr import TimeDomainVoltage
from .noiseomegaexpr import AngularFourierDomainNoiseVoltage
from .noisefexpr import FourierDomainNoiseVoltage
from .nexpr import DiscreteTimeDomainVoltage
from .kexpr import DiscreteFourierDomainVoltage
from .zexpr import ZDomainVoltage
from .phasor import PhasorDomainVoltage
from .units import u as uu


def Vname(name, kind, cache=False):

    if kind in ('s', 'laplace'):    
        return LaplaceDomainVoltage(name + '(s)')
    elif kind in ('t', 'time'):
        return TimeDomainVoltage(name.lower() + '(t)')
    elif kind in (omega0sym, omega0, 'ac'):
        return PhasorDomainVoltage(name + '(omega_0)')
    # Not caching is a hack to avoid conflicts of Vn1 with Vn1(s) etc.
    # when using subnetlists.  The alternative is a proper context
    # switch.  This would require every method to set the context.
    return voltage(name, cache=cache)            


def Vtype(kind):
    
    if isinstance(kind, str) and kind[0] == 'n':
        return AngularFourierDomainNoiseVoltage
    try:
        return {'ivp' : LaplaceDomainVoltage,
                's' : LaplaceDomainVoltage,
                'n' : AngularFourierDomainNoiseVoltage,
                'ac' : PhasorDomainVoltage,
                'dc' : ConstantVoltage,
                't' : TimeDomainVoltage,
                'time' : TimeDomainVoltage}[kind]
    except KeyError:
        return PhasorDomainVoltage


def voltage(arg, **assumptions):

    expr1 = expr(arg, **assumptions)

    try:
        expr1 = expr1.as_voltage()
    except:        
        raise ValueError('Cannot represent %s(%s) as voltage' % (expr1.__class__.__name__, expr1))
    
    return expr1.apply_unit(uu.volts)
