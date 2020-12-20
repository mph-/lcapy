class ImpedanceMixin(object):

    wrapper = 'impedance'
    is_always_causal = True

    @property
    def R(self):
        """Resistance."""
        return self.real

    @property
    def resistance(self):
        """Resistance."""
        return self.R

    @property
    def X(self):
        """Reactance."""
        return self.imag

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
    def Y(self):
        """Admittance."""

        from .admittance import admittance

        return admittance(1 / self)

    @property
    def Z(self):
        """Impedance."""
        return self
    

