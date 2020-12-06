"""This module provides the Immittance class, the base class for
Admittance and Impedance and the rest of the menagerie.

Copyright 2019--2020 Michael Hayes, UCECE

"""

from .expr import expr
from .cexpr import cExpr
from .sexpr import sExpr
from .fexpr import fExpr
from .omegaexpr import omegaExpr
from .symbols import j, omega, jomega, s, f, pi
from .sym import omegasym, fsym

class Immitance(sExpr):

    def __init__(self, val, kind=None, causal=True, positive=False, **assumptions):
        """Create an immittance (impedance/admittance).

        Note, by default, `positive` is False.  Thus if `val` is a
        string, any symbols specified in the string will be assumed to
        be complex.  For example, `Impedance('Z')`.  However, the
        symbols R and L in `Impedance('R + s * L')` will be assumed to
        be complex unless these symbols have been previously defined
        otherwise.  This may stymy some simplification."""
        
        val = expr(val, positive=positive)

        if isinstance(val, sExpr) and kind is None:
            kind = 's'
        elif isinstance(val, omegaExpr) and kind is None:
            kind = omegasym
            val = val.subs(omega, s / j)
        elif isinstance(val, fExpr) and kind is None:
            kind = fsym
            val = val.subs(f, s / (j * 2 * pi))            
        
        super(Immitance, self).__init__(val, causal=causal, **assumptions)
        self.kind = kind

    @property
    def _pexpr(self):

        """Return expression for printing."""

        kind = self.kind
        if kind is None:
            # Default to printing impedance as Z(omega)
            kind = omega
        
        return self.selectexpr(kind)
        
    def selectexpr(self, kind=None):

        if kind is None:
            kind = self.kind

        expr = self.expr
            
        if kind in ('s', 'ivp', 'super', 'laplace'):
            return expr
        elif kind in ('dc', 'time'):
            return expr.subs(self.var, 0)
        elif kind in (omegasym, omega, 'ac') or (isinstance (kind, str)
                                                 and kind.startswith('n')):
            return expr.subs(self.var, jomega.expr)
        elif kind == fsym:
            return expr.subs(self.var, (2 * pi * f * j).expr)
        return expr.subs(self.var, j * kind)

    def new(self, kind):

        return self.__class__(self.selectexpr(kind), kind)

    @property
    def R(self):
        """Resistance."""
        from .resistance import Resistance
        
        return Resistance(self.Zw.real)

    @property
    def X(self):
        """Reactance."""
        from .reactance import Reactance
        
        return Reactance(self.Zw.imag)

    @property
    def G(self):
        """Conductance."""
        from .conductance import Conductance
        
        return Conductance(self.Yw.real)

    @property
    def B(self):
        """Susceptance."""
        from .susceptance import Susceptance
        
        return Susceptance(-self.Yw.imag)   

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
    
