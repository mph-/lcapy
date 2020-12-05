"""This module provides the ImmittanceMixin class.  It provides common
methods for immitances.

Copyright 2020 Michael Hayes, UCECE

"""

from .omegaexpr import omegaExpr, Zomega, Yomega
from .sexpr import Zs, Ys
from .symbols import j, omega, jomega, s

class ImmitanceMixin(object):

    @property
    def R(self):
        """Resistance."""
        return Resistance(self.Zw.real)

    @property
    def resistance(self):
        """Resistance."""
        return self.R

    @property
    def X(self):
        """Reactance."""
        return Reactance(self.Zw.imag)

    @property
    def reactance(self):
        """Reactance."""
        return self.X

    @property
    def G(self):
        """Conductance.
        
        Note Y = G + j * B = 1 / Z = 1 / (R + j * X)
        and so G = R / (R**2 + X**2).

        Thus for DC, when X = 0, then G = 1 / R and is infinite for R
        = 0.  However, if Z is purely imaginary, i.e, R = 0 then G = 0,
        not infinity as might be expected.

        """
        return Conductance(self.Yw.real)

    @property
    def conductance(self):
        """Conductance."""
        return self.G    
    
    @property
    def B(self):
        """Susceptance."""
        return Susceptance(-self.Yw.imag)

    @property
    def susceptance(self):
        """Susceptance."""
        return self.B

    @property
    def Y(self):
        """Admittance."""
        return self.admittance

    @property
    def Z(self):
        """Impedance."""
        return self.impedance

    @property
    def Yw(self):
        """Admittance  Y(omega)."""
        return Yomega(self.admittance.selectexpr(omega))

    @property
    def Zw(self):
        """Impedance  Z(omega)."""
        return Zomega(self.impedance.selectexpr(omega))

    @property
    def Ys(self):
        """Generalized admittance  Y(s)."""
        return Ys(self.admittance.selectexpr(s))

    @property
    def Zs(self):
        """Generalized impedance  Z(s)."""
        return Zs(self.impedance.selectexpr(s))

from .resistance import Resistance
from .reactance import Reactance
from .conductance import Conductance
from .susceptance import Susceptance

