"""This file provides the immittance class for admittances and
impedances and the rest of the menagerie.

Copyright 2019 Michael Hayes, UCECE

"""

from .cexpr import cExpr
from .sexpr import sExpr
from .omegaexpr import omegaExpr
from .symbols import j, omega, jomega, s
from .sym import omegasym

class Immitance(sExpr):
    
    def __init__(self, val, kind=None, causal=True, **assumptions):

        super(Immitance, self).__init__(val, causal=causal, **assumptions)
        self.kind = kind

    @property
    def pexpr(self):

        kind = self.kind
        if kind is None:
            # Default to printing impedance as Z(omega)
            kind = omega
        
        return self.select(kind).expr
        
    @property            
    def jomega(self):
        return self(jomega)

    @property    
    def s(self):
        return self(s)    

    def select(self, kind=None):

        if kind is None:
            kind = self.kind

        if kind in ('s', 'ivp', 'super', 'laplace'):
            return self
        elif kind in ('dc', 'time'):
            return cExpr(self.subs(0))
        elif isinstance(kind, str) and kind[0] == 'n':
            return self(jomega)
        elif kind in (omegasym, omega, 'ac'):
            return self(jomega)
        return omegaExpr(self.subs(j * kind))    

    @property
    def R(self):
        """Resistance."""
        return self.Zw.real

    @property
    def X(self):
        """Reactance."""
        return self.Zw.imag

    @property
    def G(self):
        """Conductance."""
        return self.Yw.real

    @property
    def B(self):
        """Susceptance."""
        return -self.Yw.imag    

    
class ImmitanceMixin(object):

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
