"""This module provides the ImmittanceMixin class.  It provides common
methods for immittances.

Copyright 2021 Michael Hayes, UCECE

"""

class ImmittanceMixin(object):

    @property
    def R(self):
        """Resistance."""
        from .symbols import s
        
        ret = self.Z.real
        if self.is_causal:
            ret = ret.replace(s.real, 0)
        return ret

    @property
    def resistance(self):
        """Resistance."""
        return self.R

    @property
    def is_lossy(self):
        """Has some none zero resistance component."""
        return self.R != 0

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
        """Conductance.
        
        Note Y = G + j * B = 1 / Z = 1 / (R + j * X)
        and so G = R / (R**2 + X**2).

        Thus for DC, when X = 0, then G = 1 / R and is infinite for R
        = 0.  However, if Z is purely imaginary, i.e, R = 0 then G = 0,
        not infinity as might be expected.

        """
        from .symbols import s
        
        ret = self.Y.real
        if self.is_causal:
            ret = ret.replace(s.real, 0)
        return ret        

    @property
    def conductance(self):
        """Conductance."""
        return self.G    
    
    @property
    def B(self):
        """Susceptance."""
        return -self.Y.imag

    @property
    def Y(self):
        """Admittance.

        The admittance is expressed in jomega form for AC circuits
        and in s-domain for for transient circuits.
        
        Use Y(omega) or Y(s) to achieve the desired form."""
        
        return self.admittance

    @property
    def Z(self):
        """Impedance.

        The impedance is expressed in jomega form for AC circuits
        and in s-domain for for transient circuits.  

        Use Z(omega) or Z(s) to achieve the desired form."""

        return self.impedance    

    def network(self, form='default'):
        """Synthesise a network with an equivalent impedance.
        `form` includes: cauerI, cauerII, fosterI, fosterII.

        Note some methods generate networks with negative value
        components."""

        from .synthesis import network
        
        return network(self.Z, form)
