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


Copyright 2014--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .acdc import ACChecker
from .domains import PhasorTimeDomain, PhasorFrequencyDomain
from .sym import j, omegasym
from .expr import expr
from .functions import sin, cos, exp, sqrt
from .expr import Expr
from .omegaexpr import omega
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



class PhasorExpression(Expr):

    def __init__(self, val, **assumptions):

        if isinstance(val, PhasorExpression):
            assumptions['omega'] = val.omega
        elif 'omega' not in assumptions or assumptions['omega'] is None:
            assumptions['omega'] = omegasym
        
        assumptions['ac'] = True
        super (PhasorExpression, self).__init__(val, **assumptions)
    
    @property
    def omega(self):
        """Return angular frequency."""

        return self.assumptions.get('omega', None)

    @property
    def var(self):
        """Return angular frequency."""

        return self.omega

    def fourier(self, **assumptions):
        """Fourier transform."""

        return self.time().fourier()

    def angular_fourier(self, **assumptions):
        """Angular Fourier transform."""

        if self.has(omega):
            print('Warning: expression contains omega, should substitute with a different symbol')
        
        return self.time().angular_fourier()            
    
    def laplace(self, **assumptions):
        """Convert to Laplace domain representation."""

        return self.time().laplace()

    def phasor(self):
        """Convert to phasor representation."""
        
        return self.__class__(self, **self.assumptions)

    def rms(self):
        """Return root mean square."""
        
        return abs(self) * sqrt(2) / 2

    def plot(self, wvector=None, **kwargs):
        """Plot polar diagram for a time-domain phasor or frequency response
        for frequency-domain phasor.  For the latter, wvector
        specifies the angular frequencies.  If it is a tuple, it sets
        the angular frequency limits."""
        
        from .plot import plot_phasor, plot_angular_frequency

        if self.is_phasor_frequency_domain:
            return plot_angular_frequency(self, wvector, **kwargs)            
        
        return plot_phasor(self, **kwargs)

    def _mul_compatible(self, x):

        if not hasattr(x, 'omega'):
            return True

        if self.omega == x.omega:
            return True
        
        raise ValueError('Incompatible phasor angular frequencies %s and %s' %
                         (self.omega, x.omega))
        

class PhasorTimeDomainExpression(PhasorTimeDomain, PhasorExpression):
    """This is a phasor domain base class for voltages and currents."""

    is_phasor_domain = True
    is_phasor_time_domain = True

    def as_expr(self):
        return PhasorTimeDomainExpression(self)

    @classmethod
    def from_time(cls, expr, omega=None, **assumptions):

        from .symbols import t

        if expr.is_admittance or expr.is_impedance or expr.is_transfer:
            print('Should convert %s expression to Laplace-domain first' % expr.quantity)

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

        if check.omega == 0:
            print('Warning, DC phasor.')
        
        result = check.amp * exp(j * check.phase)
        assumptions['omega'] = check.omega

        return cls.change(expr, result, domain='phasor', **assumptions)

    def time(self, **assumptions):
        """Convert to time domain representation."""
        from .symbols import t
        
        omega1 = self.omega
        if isinstance(omega1, Expr):
            # TODO: Fix inconsistency.  Sometimes omega is a symbol.
            omega1 = omega1.expr

        if self.is_complex:
            result = self.expr * exp(j * omega1 * t)
        else:
            result = self.real.expr * cos(omega1 * t) - self.imag.expr * sin(omega1 * t)

        return TimeDomainExpression(result).as_quantity(self.quantity)
    

class PhasorFrequencyDomainExpression(PhasorFrequencyDomain, PhasorExpression):
    """This represents the ratio of two-phasors; for example
    an impedance, an admittance, or a transfer function."""

    is_phasor_frequency_domain = True    
    is_transform_domain = True

    @classmethod
    def from_laplace(cls, expr, omega=None, **assumptions):

        from .symbols import s        

        if omega is None:
            omega = omegasym

        if expr.is_voltage or expr.is_current:
            print('Should convert %s expression to time-domain first' % expr.quantity)

        result = expr.laplace(**assumptions).replace(s, j * omega)
        return cls.change(expr, PhasorFrequencyDomainExpression(result, omega=omega,
                                                                **assumptions))
    
    def time(self, **assumptions):
        """Convert to time domain representation."""
        from .symbols import s
        from .sexpr import LaplaceDomainExpression

        omega = self.omega
        result = LaplaceDomainExpression(self.replace(j * omega, s))
        return self.change(result.time())

    def as_expr(self):
        return PhasorFrequencyDomainExpression(self)
    

def phasor(arg, omega=None, **assumptions):
    """Create phasor.   

    If arg has the form A * cos(w * t + phi) the phasor 
    A * exp(j * phi) of angular frequency w is returned.

    If arg is a constant C, the phasor C is created with omega as the
    angular frequency (this defaults to 'omega' if not specified).

    """

    arg = expr(arg)

    if arg.is_time_domain:
        # Expecting AC signal.
        return PhasorTimeDomainExpression.from_time(arg, omega=omega, **assumptions)
    elif arg.is_unchanging:
        # Expecting phasor (complex amplitude)        
        return PhasorTimeDomainExpression(arg, omega=omega, **assumptions)
    else:
        # Is this sensible?
        return PhasorFrequencyDomainExpression(arg, omega=omega, **assumptions)        
    

from .expressionclasses import expressionclasses

expressionclasses.register('phasor', PhasorTimeDomainExpression,
                           PhasorFrequencyDomainExpression,
                           ('voltage', 'current', 'voltagesquared', 'currentsquared'))
from .texpr import TimeDomainExpression
from .expr import Expr

