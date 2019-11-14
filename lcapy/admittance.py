"""This file provides the Admittance class.  This is a generalized
impedance (s-domain) and converts to other representations.

Copyright 2019 Michael Hayes, UCECE

"""

from .symbols import s
from .immitance import Immitance

class Admittance(Immitance):    
    """Generic admittance class."""

    @property
    def Yw(self):
        return self.jomega

    @property
    def Zw(self):
        return 1 / self.Yw

    def cpt(self):
        """Create oneport component."""

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
