"""This file provides the immittance class for admittances and
impedances and the rest of the menagerie.

Copyright 2019 Michael Hayes, UCECE

"""

class Immitance(object):

    @property
    def R(self):
        """Resistance."""
        return self.Z.real

    @property
    def resistance(self):
        """Resistance."""
        return self.R

    @property
    def X(self):
        """Reactance."""
        return self.Z.imag

    @property
    def reactance(self):
        """Reactance."""
        return self.X    

    @property
    def G(self):
        """Conductance."""
        return self.Y.real

    @property
    def conductance(self):
        """Conductance."""
        return self.G    
    
    @property
    def B(self):
        """Susceptance."""
        return -self.Y.imag

    @property
    def susceptance(self):
        """Susceptance."""
        return self.B

    @property
    def admittance(self):
        """Admittance."""
        return self.Y

    @property
    def impedance(self):
        """Impedance."""
        return self.Z

    @property
    def generalized_admittance(self):
        """Generalized admittance."""
        return self.Ys

    @property
    def generalized_impedance(self):
        """Generalized impedance."""
        return self.Zs
    
