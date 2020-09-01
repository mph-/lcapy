"""This module provides the Admittance class.  This is a generalized
impedance (s-domain) and converts to other representations.

Copyright 2019-2020 Michael Hayes, UCECE

"""
from __future__ import division
from .symbols import s
from .immitance import Immitance

class Admittance(Immitance):    
    """Generic admittance class.

    Y(omega) = G(omega) + j * B(omega)

    where G is the conductance and B is the susceptance.

    Admittance is the reciprocal of impedance,

    Z(omega) = 1 / Y(omega)

    """

    def __rtruediv__(self, x):
        """Reverse true divide"""

        # TODO: handle Current / Admittance -> Voltage etc.
        from .impedance import Impedance
        return Impedance(x / self.expr)

    @property
    def Y(self):
        return self    
    
    @property
    def Z(self):
        return 1 / self

    @property
    def Yw(self):
        return self.jomega

    @property
    def Zw(self):
        return 1 / self.Yw

    def cpt(self):
        """Create oneport component.  See also network."""        

        from .oneport import G, C, L, Y

        if self.is_number or self.is_dc:
            return G(self.expr)

        y = self * s

        if y.is_number:
            return L((1 / y).expr)

        y = self / s

        if y.is_number:
            return C(y.expr)

        return Y(self)    

    def network(self, form='default'):
        """Synthesise a network with an equivalent admittance.
        `form` includes: cauerI, cauerII, fosterI, fosterII.

        Note some methods generate networks with negative value
        components."""
        
        return self(s).network(form)
    
    
