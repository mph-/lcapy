"""This module provides the PhasorDomain classes to represent phasors
 for AC analysis.

A phasor represents the amplitude and phase for a single sinusoid.  By
default the angular frequency is omega_0 but it can be any number or
symbol.

Phasors straddle the time and frequency domains and this is my excuse
for the confusing phasor classes.

A phasor is described by an amplitude and phase.  There is also an
implicit frequency.  The amplitude and phase are usually functions of
frequency.  A phasor is useful to describe an AC voltage or current
signal.

A ratio of two phasors (of the same frequency) is no longer a
time-domain signal but a frequency domain quantity.  The frequency
dependence is explicit.  A phasor ratio is useful to describe an
impedance, immitance, or transfer function.  These are frequency
domain concepts.  A phasor ratio can be inferred from the Laplace
domain by substituting jomega for s, where omega is the angular
frequency of the phasor.


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

__all__ = ('phasor', )

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

    is_phasor_domain = True
    domain_name = 'Phasor'

    def __init__(self, val, **assumptions):

        if isinstance(val, PhasorDomainExpression):
            assumptions['omega'] = val.omega
        elif 'omega' not in assumptions:
            assumptions['omega'] = omega0sym                    

        assumptions['ac'] = True
        super (PhasorDomainExpression, self).__init__(val, **assumptions)
    
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

        return self.assumptions.get('omega', None)

    @property
    def var(self):
        """Return angular frequency."""

        return self.omega

    def __compat_mul__(self, x, op):

        cls = self.__class__
        xcls = x.__class__

        # Perhaps check explicitly for int, float?
        if not isinstance(x, Expr):
            return cls, self, cls(x), self.assumptions

        if isinstance(x, (AngularFourierDomainExpression, ConstantExpression)):
            return cls, self, x, self.assumptions

        if not isinstance(x, PhasorDomainTimeExpression):
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


class PhasorDomainTimeExpression(PhasorDomainExpression):
    """This is a phasor domain base class for voltages and currents."""

    is_phasor_domain = True
    is_phasor_time_domain = True

    def as_expr(self):
        return PhasorDomainTimeExpression(self)

    @classmethod
    def from_time(cls, expr, omega=None, **assumptions):

        from .symbols import t

        assumptions['ac'] = True

        if expr.is_transform_domain:
            print('Warning, converting %s-domain to time-domain first.' %
                  expr.domain)
            expr = expr.time()
            
        check = ACChecker(expr, t)
        if not check.is_ac:
            raise ValueError(
                'Do not know how to convert %s to phasor.  Expecting an AC signal.' % expr)
        
        if omega is not None and check.omega != omega.expr:
            raise ValueError('Expecting omega=%s, found omega=%s.' % (omega, check.omega))
        result = check.amp * exp(j * check.phase)
        assumptions['omega'] = check.omega

        return cls.wrap(expr, PhasorDomainTimeExpression(result, **assumptions))

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

        return TimeDomainExpression(result)
    

class PhasorDomainFrequencyExpression(PhasorDomainExpression):
    """This represents the ratio of two-phasors; for example
    an impedance, an admittance, or a transfer function."""

    is_phasor_frequency_domain = True    
    is_transform_domain = True

    @classmethod
    def from_laplace(cls, expr, omega=None, **assumptions):

        from .symbols import s        

        if omega is None:
            raise ValueError('omega unspecified for conversion of %s-domain to phasor-domain' % expr.domain)

        result = expr.laplace(**assumptions).replace(s, j * omega)
        return cls.wrap(expr, PhasorDomainFrequencyExpression(result, omega=omega,
                                                              **assumptions))
    
    def time(self, **assumptions):
        """Convert to time domain representation."""
        from .symbols import jw, s

        return self.wrap(cls, TimeDomainExpression(self.replace(jw, s)).time(causal=True))

    def as_expr(self):
        return PhasorDomainFrequencyExpression(self)
    

class PhasorDomainAdmittance(AdmittanceMixin, PhasorDomainFrequencyExpression):
    """PhasorDomain admittance"""
    pass


class PhasorDomainImpedance(ImpedanceMixin, PhasorDomainFrequencyExpression):
    """PhasorDomain impedance"""
    pass


class PhasorDomainTransferFunction(TransferMixin, PhasorDomainFrequencyExpression):
    """PhasorDomain transfer function response."""
    pass

    
class PhasorDomainVoltage(VoltageMixin, PhasorDomainTimeExpression):
    """t-domain voltage (units V) parameterized as a phasor
    of a single angular frequency, omega0."""
        
    def cpt(self):
        from .oneport import Vac
        return Vac(self, 0, self.omega)            

    
class PhasorDomainCurrent(CurrentMixin, PhasorDomainTimeExpression):
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
             (PhasorDomainTimeExpression, PhasorDomainAdmittance): (None, PhasorDomainImpedance),
             (PhasorDomainTimeExpression, PhasorDomainImpedance): (None, PhasorDomainAdmittance)}

    
def phasor(arg, **assumptions):
    """Create phasor."""

    arg = expr(arg)

    return PhasorDomainTimeExpression(arg, **assumptions)
    
    
from .sexpr import LaplaceDomainExpression
from .texpr import TimeDomainExpression, TimeDomainVoltage, TimeDomainCurrent
from .expr import Expr

