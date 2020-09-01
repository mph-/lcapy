"""This module provides the Impedance class.  This is a generalized
impedance (s-domain) and converts to other representations.

Copyright 2019-2020 Michael Hayes, UCECE

"""
from __future__ import division
from .symbols import s
from .immitance import Immitance

class Impedance(Immitance):
    """Generic impedance class.

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

    def cpt(self):
        """Create oneport component.  See also network."""
        
        from .oneport import R, C, L, Z

        if self.is_number or self.is_dc:
            return R(self.expr)

        z = self * s

        if z.is_number:
            return C((1 / z).expr)

        z = self / s

        if z.is_number:
            return L(z.expr)

        return Z(self)

    def network(self, form='default'):
        """Synthesise a network with an equivalent impedance.
        `form` includes: cauerI, cauerII, fosterI, fosterII.

        Note some methods generate networks with negative value
        components."""
        
        return self(s).network(form)
