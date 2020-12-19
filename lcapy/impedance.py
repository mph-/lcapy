"""This module provides the Impedance class.  This is a generalized
impedance (s-domain) and converts to other representations.

Copyright 2019--2020 Michael Hayes, UCECE

"""
from __future__ import division
from .expr import expr
from .fexpr import FourierDomainExpression, FourierDomainImpedance
from .omegaexpr import AngularFourierDomainExpression, AngularFourierDomainImpedance
from .sexpr import LaplaceDomainExpression, LaplaceDomainImpedance
from .symbols import s
from .immittance import Immittance
from .units import u as uu


class Impedance(Immittance):
    """Generic impedance class.

    This represents both the phasor impedance, Z(omega), and the
    s-domain impedance (sometimes called generalized imepdance), Z(s).
    Unfortunately, both are called impedance although Z(omega) is more
    common.

    By default the representation uses the form that the impedance was
    defined.  Otherwise, use Z(s) or Z(omega) to get the desired form.

    Z(omega) = R(omega) + j * X(omega)

    where R is the resistance and X is the reactance.

    Impedance is the reciprocal of admittance,

    Y(omega) = 1 / Z(omega)

    """

    def __rtruediv__(self, x):
        """Reverse true divide"""

        # TODO: handle Voltage / Impedance -> Current etc.
        from .admittance import Admittance
        return Admittance(x / self.expr)

    @property
    def Y(self):
        return 1 / self
    
    @property
    def Z(self):
        return self    
    
    @property
    def Yw(self):
        return 1 / self.Zw

    @property
    def Zw(self):
        return self.jomega

    @property
    def Ys(self):
        return 1 / self.Zs

    @property
    def Zs(self):
        return self(s)
    
    def cpt(self):
        """Create oneport component.  See also network."""
        
        from .oneport import R, C, L, Z

        if self.is_number or self.is_dc:
            return R(self.expr)

        z = self.Zs * s

        if z.is_number:
            return C((1 / z).expr)

        z = self.Zs / s

        if z.is_number:
            return L(z.expr)

        return Z(self.orig)

    def network(self, form='default'):
        """Synthesise a network with an equivalent impedance.
        `form` includes: cauerI, cauerII, fosterI, fosterII.

        Note some methods generate networks with negative value
        components."""
        
        return self(s).network(form)


def impedance(arg):

    mapping = {ConstantExpression: AngularFourierDomainImpedance,
               TimeDomainExpression: TimeDomainImpedance,
               LaplaceDomainExpression: LaplaceDomainImpedance,
               FourierDomainExpression: FourierDomainImpedance,
               AngularFourierDomainExpression: AngularFourierDomainImpedance}
    
    expr1 = expr(arg)
    if expr1.__class__ in mapping:
        expr1 = mapping(expr1)
    
    return expr1.apply_unit(uu.ohms)    
    
