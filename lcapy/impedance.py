"""This module provides the Impedance class.  This is a generalized
impedance (s-domain) and converts to other representations.

Copyright 2019 Michael Hayes, UCECE

"""


from .symbols import s
from .immitance import Immitance

class Impedance(Immitance):
    """Generic impedance class.

    Z(omega) = R(omega) + j * X(omega)

    where R is the resistance and X is the reactance.

    Impedance is the reciprocal of admittance,

    Y(omega) = 1 / Z(omega)

    """

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
        """Create oneport component."""
        
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
    

