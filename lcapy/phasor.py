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
from .domains import PhasorDomain
from .sym import j, omegasym
from .expr import expr
from .functions import sin, cos, exp, sqrt
from .expr import Expr
from .cexpr import ConstantExpression
from .omegaexpr import AngularFourierDomainExpression, omega
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



class PhasorDomainExpression(PhasorDomain, Expr):

    def __init__(self, val, **assumptions):

        if isinstance(val, PhasorDomainExpression):
            assumptions['omega'] = val.omega
        elif 'omega' not in assumptions:
            assumptions['omega'] = omegasym

        assumptions['ac'] = True
        super (PhasorDomainExpression, self).__init__(val, **assumptions)
    
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

    def _mul_compatible(self, x):

        if not hasattr(x, 'omega'):
            return True

        if self.omega == x.omega:
            return True
        
        raise ValueError('Incompatible phasor angular frequencies %s and %s' %
                         (self.omega, x.omega))
        

class PhasorDomainTimeExpression(PhasorDomainExpression):
    """This is a phasor domain base class for voltages and currents."""

    is_phasor_domain = True
    is_phasor_time_domain = True

    def as_expr(self):
        return PhasorDomainTimeExpression(self)

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
            omega = omegasym

        if expr.is_voltage or expr.is_current:
            print('Should convert %s expression to time-domain first' % expr.quantity)

        result = expr.laplace(**assumptions).replace(s, j * omega)
        return cls.change(expr, PhasorDomainFrequencyExpression(result, omega=omega,
                                                                **assumptions))
    
    def time(self, **assumptions):
        """Convert to time domain representation."""
        from .symbols import s

        omega = self.omega
        result = LaplaceDomainExpression(self.replace(j * omega, s))
        return self.change(result.time())

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
    """Phasor-domain voltage (units V) parameterized as a phasor
    of a single angular frequency."""
        
    def cpt(self):
        from .oneport import Vac
        return Vac(self, 0, self.omega)            

    
class PhasorDomainCurrent(CurrentMixin, PhasorDomainTimeExpression):
    """Phasor-domain current (units V) parameterized as a phasor
    of a single angular frequency."""    

    def cpt(self):
        from .oneport import Iac
        return Iac(self, 0, self.omega)


def phasor(arg, **assumptions):
    """Create phasor."""

    arg = expr(arg)

    if arg.is_time_domain:
        return PhasorDomainTimeExpression(arg, **assumptions)
    else:
        return PhasorDomainFrequencyExpression(arg, **assumptions)
    

from .expressionclasses import expressionclasses

classes = expressionclasses.make(PhasorDomainExpression)

classes['voltage'] = PhasorDomainVoltage
classes['current'] = PhasorDomainCurrent
classes['admittance'] = PhasorDomainAdmittance
classes['impedance'] = PhasorDomainImpedance
classes['transfer'] = PhasorDomainTransferFunction
expressionclasses.add('phasor', classes)

from .sexpr import LaplaceDomainExpression
from .texpr import TimeDomainExpression, TimeDomainVoltage, TimeDomainCurrent
from .expr import Expr

