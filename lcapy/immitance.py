"""This file provides the Immittance class, the base class for
Admittance and Impedance and the rest of the menagerie.

Copyright 2019 Michael Hayes, UCECE

"""

from .expr import expr
from .cexpr import cExpr
from .sexpr import sExpr
from .omegaexpr import omegaExpr
from .symbols import j, omega, jomega, s
from .sym import omegasym

class Immitance(sExpr):
    
    def __init__(self, val, kind=None, causal=True, **assumptions):

        val = expr(val)
        if isinstance(val, omegaExpr) and kind is None:
            val = val.subs(omega, s / j)
        
        super(Immitance, self).__init__(val, causal=causal, **assumptions)
        self.kind = kind

    @property
    def _pexpr(self):
        """Return expression for printing."""

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

    def new(self, kind):

        return self.__class__(self.select(kind), kind)

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

    def oneport(self):
        """Create oneport component."""
        return self.cpt()

    @property
    def real(self):
        """Real part of Z(omega)."""
        return self.Zw.real

    @property
    def imag(self):
        """Imaginary part of Z(omega)."""
        return self.Zw.imag

    @property
    def abs(self):
        """Absolute part of Z(omega)."""
        return self.Zw.abs

    @property
    def phase(self):
        """Phase of Z(omega) (radians)."""
        return self.Zw.phase

    @property
    def phase_degrees(self):
        """Phase of Z(omega) (degrees)."""
        return self.Zw.phase_degrees        
    
    
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
        """Conductance.
        
        Note Y = G + j * B = 1 / Z = 1 / (R + j * X)
        and so G = R / (R**2 + X**2).

        Thus for DC, when X = 0, then G = 1 / R and is infinite for R
        = 0.  However, if Z is purely imaginary, i.e, R = 0 then G = 0,
        not infinity as might be expected.

        """
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
        return self.admittance.jomega

    @property
    def Zw(self):
        """Impedance  Z(omega)."""
        return self.impedance.jomega        

    @property
    def Ys(self):
        """Generalized admittance  Y(s)."""
        return self.admittance.s

    @property
    def Zs(self):
        """Generalized impedance  Z(s)."""
        return self.impedance.s
