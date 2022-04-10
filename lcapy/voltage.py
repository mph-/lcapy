"""This module provides voltage support.

Copyright 2020 Michael Hayes, UCECE

"""
from .expr import expr
from .classmap import domain_kind_to_symbol, domain_kind_quantity_to_class


def Vname(name, kind, cache=False):

    # Not caching is a hack to avoid conflicts of Vn1 with Vn1(s) etc.
    # when using subnetlists.  The alternative is a proper context
    # switch.  This would require every method to set the context.

    if kind == 't':
        name = name.lower()

    undef = domain_kind_to_symbol(kind, name)
    cls = domain_kind_quantity_to_class(kind, 'voltage')
    return cls(undef, cache=cache)


def Vtype(kind):

    return domain_kind_quantity_to_class(kind, 'voltage')


def voltage(arg, **assumptions):

    expr1 = expr(arg, **assumptions)

    if 'nid' in assumptions:
        from .noisefexpr import FourierNoiseDomainVoltage
        from .noiseomegaexpr import AngularFourierNoiseDomainVoltage

        if expr1.is_fourier_domain or expr1.is_constant_domain:
            expr1 = FourierNoiseDomainVoltage(expr1)
        elif expr1.is_angular_fourier_domain:
            expr1 = AngularFourierNoiseDomainVoltage(expr1)
        else:
            raise ValueError(
                'Cannot represent noise voltage in %s domain' % expr1.domain)

    try:
        expr1 = expr1.as_voltage()
    except:
        raise ValueError('Cannot represent %s(%s) as voltage' %
                         (expr1.__class__.__name__, expr1))

    return expr1


def noisevoltage(arg, **assumptions):
    """Create a new noise voltage with specified amplitude spectral density."""

    nid = assumptions.get('nid', None)
    positive = assumptions.get('positive', True)

    return voltage(arg, nid=None, positive=positive, **assumptions)


def phasorvoltage(arg, omega=None, **assumptions):
    from .phasor import phasor

    return phasor(arg, omega, **assumptions).as_voltage()
