"""This module provides the PhasorDomainExpression class to represent phasors
 for AC analysis.

A phasor represents the amplitude and phase for a single sinusoid.  By
default the angular frequency is omega_0 but it can be any number or
symbol.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .acdc import ACChecker
from .sym import j, omega0sym
from .expr import expr
from .functions import sin, cos, exp, sqrt
from .expr import Expr
from .cexpr import ConstantExpression
from .omegaexpr import AngularFourierDomainExpression
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
from .admittancemixin import AdmittanceMixin
from .impedancemixin import ImpedanceMixin
from .transfermixin import TransferMixin


# The phasor domain is different from the Fourier and Laplace domain
# since there is an implicit anglular frequency.  This is only needed
# for voltage and current (vorrent) expressions.

# The phasor domain immittance can be found from the Laplace domain
# impedance/admittance by substituting s with j * omega.  In many
# cases the result is the same as the Fourier domain immittance but
# not always.  The most important case is the impedance of a
# capacitor.  In the Laplace domain it is 1/(s C); in the phasor
# domain it is 1/(j omega C); in the Fourier domain it is 1/(2 j omega
# C) + delta(omega) / (2 C).



class PhasorDomainExpression(Expr):

    is_phasor = True
    domain_name = 'Phasor'
    
    @classmethod
    def make(cls, expr, **assumptions):
        from .symbols import s, jw, t

        if expr.is_phasor:
            return expr.wrap(expr)
        
        if expr.is_voltage or expr.is_current:
        
            check = ACChecker(expr.time(), t)
            if not check.is_ac:
                raise ValueError('Do not know how to convert %s to phasor' % expr)
            phasor = PhasorDomainVorrent(check.amp * exp(j * check.phase), omega=check.omega)
            return expr.wrap(phasor)
        else:
            result = expr.wrap(PhasorDomainRatio(expr.laplace().replace(s, jw)))
        return result

    def _class_by_quantity(self, quantity):

        if quantity == 'voltage':
            return PhasorDomainVoltage
        elif quantity == 'current':
            return PhasorDomainCurrent
        elif quantity == 'impedance':
            return PhasorDomainImpedance
        elif quantity == 'admittance':
            return PhasorDomainAdmittance
        elif quantity == 'transfer':
            return PhasorDomainTransferFunction                
        raise ValueError('Unknown quantity %s' % quantity)
    
    def as_expr(self):
        return PhasorDomainExpression(self)

    def as_voltage(self):
        return PhasorDomainVoltage(self)

    def as_current(self):
        return PhasorDomainCurrent(self)    

    def as_impedance(self):
        return PhasorDomainImpedance(self)

    def as_admittance(self):
        return PhasorDomainAdmittance(self)

    def as_transfer(self):
        return PhasorDomainTransferFunction(self)    
    
    @property
    def omega(self):
        """Return angular frequency."""

        return self.assumptions['omega']

    @property
    def var(self):
        """Return angular frequency."""

        return self.omega

    def __compat_add__(self, x, op):

        cls = self.__class__
        xcls = x.__class__

        # Special case for zero.
        if isinstance(x, int) and x == 0:
            return cls, self, cls(x), self.assumptions

        if not isinstance(x, PhasorDomainExpression):
            raise TypeError('Incompatible arguments %s and %s for %s' %
                            (repr(self), repr(x), op))

        if self.omega != x.omega:
            raise ValueError('Cannot combine %s(%s, omega=%s)'
                             ' with %s(%s, omega=%s)' %
                             (cls.__name__, self, self.omega,
                              xcls.__name__, x, x.omega))
        return cls, self, x, self.assumptions

    def __compat_mul__(self, x, op):

        cls = self.__class__
        xcls = x.__class__

        # Perhaps check explicitly for int, float?
        if not isinstance(x, Expr):
            return cls, self, cls(x), self.assumptions

        if isinstance(x, (AngularFourierDomainExpression, ConstantExpression)):
            return cls, self, x, self.assumptions

        if not isinstance(x, PhasorDomainExpression):
            raise TypeError('Incompatible arguments %s and %s for %s' %
                            (repr(self), repr(x), op))

        if self.omega != x.omega:
            raise ValueError('Cannot combine %s(%s, omega=%s)'
                             ' with %s(%s, omega=%s)' %
                             (cls.__name__, self, self.omega,
                              xcls.__name__, x, x.omega))
        return cls, self, x, self.assumptions

    def fourier(self, **assumptions):
        """Fourier transform."""

        return self.time().fourier()

    def angular_fourier(self, **assumptions):
        """Angular Fourier transform."""

        return self.time().angular_fourier()            
    
    def laplace(self, **assumptions):
        """Convert to Laplace domain representation."""

        return self.time().laplace()

    @property
    def abs(self):
        """Return magnitude"""

        return expr(self.expr).abs

    @property
    def magnitude(self):
        """Return magnitude"""

        return expr(self.expr).magnitude    

    @property
    def phase(self):
        """Return phase in radians."""

        return expr(self.expr).phase

    @property
    def sign(self):
        """Return sign."""

        return expr(self.expr).sign

    def phasor(self):
        """Convert to phasor representation."""
        return self.__class__(self, **self.assumptions)

    def rms(self):
        return {PhasorDomainVoltage: TimeDomainVoltage,
                PhasorDomainCurrent : TimeDomainCurrent}[self.__class__](0.5 * self)

    def plot(self, **kwargs):

        from .plot import plot_phasor
        return plot_phasor(self, **kwargs)


class PhasorDomainVorrent(PhasorDomainExpression):
    """This is a phasor domain base class for voltages and currents."""
    
    def __init__(self, val, **assumptions):

        assumptions['ac'] = True

        if 'omega' not in assumptions:
            if hasattr(val, 'omega'):
                assumptions['omega'] = val.omega
            else:
                assumptions['omega'] = omega0sym        
        
        super (PhasorDomainExpression, self).__init__(val, **assumptions)

    def time(self, **assumptions):
        """Convert to time domain representation."""
        from .symbols import t
        
        omega = self.omega
        if isinstance(omega, Expr):
            # TODO: Fix inconsistency.  Sometimes omega is a symbol.
            omega = omega.expr
            
        if self.is_complex:
            result = self.expr * exp(j * omega * t)
        else:
            result = self.real.expr * cos(omega * t) - self.imag.expr * sin(omega * t)

        return self.wrap(TimeDomainExpression(result))


class PhasorDomainRatio(PhasorDomainExpression):
    
    def time(self, **assumptions):
        """Convert to time domain representation."""
        from .symbols import jw, s

        return self.wrap(LaplaceDomainExpression(self.replace(jw, s)).time(causal=True))
        

class PhasorDomainAdmittance(AdmittanceMixin, PhasorDomainRatio):
    """PhasorDomain admittance"""
    pass


class PhasorDomainImpedance(ImpedanceMixin, PhasorDomainRatio):
    """PhasorDomain impedance"""
    pass


class PhasorDomainTransferFunction(TransferMixin, PhasorDomainRatio):
    """PhasorDomain transfer function response."""
    pass

    
class PhasorDomainVoltage(VoltageMixin, PhasorDomainVorrent):
    """t-domain voltage (units V) parameterized as a phasor
    of a single angular frequency, omega0."""
        
    def cpt(self):
        from .oneport import Vac
        return Vac(self, 0, self.omega)            

    
class PhasorDomainCurrent(CurrentMixin, PhasorDomainVorrent):
    """t-domain current (units V) parameterized as a phasor
    of a single angular frequency, omega0."""    

    def cpt(self):
        from .oneport import Iac
        return Iac(self, 0, self.omega)


# TODO, allow PhasorDomainVoltage * PhasorDomainVoltage etc
mul_table = {(PhasorDomainVoltage, PhasorDomainAdmittance): (None, PhasorDomainCurrent),
             (PhasorDomainCurrent, PhasorDomainImpedance): (None, PhasorDomainVoltage),
             (PhasorDomainCurrent, PhasorDomainTransferFunction): (None, PhasorDomainCurrent),
             (PhasorDomainVoltage, PhasorDomainTransferFunction): (None, PhasorDomainVoltage)}

# If have a constant, etc, pass through expr first.  

div_table = {(PhasorDomainVoltage, PhasorDomainImpedance): (None, PhasorDomainCurrent),
             (PhasorDomainCurrent, PhasorDomainAdmittance): (None, PhasorDomainVoltage),
             (PhasorDomainCurrent, PhasorDomainCurrent): (None, PhasorDomainTransferFunction),
             (PhasorDomainVoltage, PhasorDomainVoltage): (None, PhasorDomainTransferFunction),
             (PhasorDomainExpression, PhasorDomainAdmittance): (None, PhasorDomainImpedance),
             (PhasorDomainExpression, PhasorDomainImpedance): (None, PhasorDomainAdmittance)}

    
def phasor(arg, **assumptions):
    """Create phasor."""

    return PhasorDomainVorrent(arg, **assumptions)
    
    
from .sexpr import LaplaceDomainExpression
from .texpr import TimeDomainExpression, TimeDomainVoltage, TimeDomainCurrent
from .expr import Expr
from .phasor import PhasorDomainExpression
