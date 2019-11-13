"""This file provides the immittance class for admittances and
impedances and the rest of the menagerie.

Copyright 2019 Michael Hayes, UCECE

"""

class Immitance(object):

    @property
    def R(self):
        """Resistance."""
        return self.Zw.real

    @property
    def resistance(self):
        """Resistance."""
        return self.R

    @property
    def X(self):
        """Reactance."""
        return self.Zw.imag

    @property
    def reactance(self):
        """Reactance."""
        return self.X    

    @property
    def G(self):
        """Conductance."""
        return self.Yw.real

    @property
    def conductance(self):
        """Conductance."""
        return self.G    
    
    @property
    def B(self):
        """Susceptance."""
        return -self.Yw.imag

    @property
    def susceptance(self):
        """Susceptance."""
        return self.B

    @property
    def admittance(self):
        """Admittance  Y(omega)."""
        return self.YY

    @property
    def impedance(self):
        """Impedance  Z(omega)."""
        return self.ZZ

    @property
    def Y(self):
        """Admittance."""
        return self.YY

    @property
    def Z(self):
        """Impedance."""
        return self.ZZ

    @property
    def Yw(self):
        """Admittance  Y(omega)."""
        return self.YY.jomega

    @property
    def Zw(self):
        """Impedance  Z(omega)."""
        return self.ZZ.jomega        

    @property
    def Ys(self):
        """Generalized admittance  Y(s)."""
        return self.YY.s

    @property
    def Zs(self):
        """Generalized impedance  Z(s)."""
        return self.ZZ.s
