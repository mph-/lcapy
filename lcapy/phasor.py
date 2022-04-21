"""This module provides the PhasorDomain classes to represent phasors
 for AC analysis.

A phasor represents the amplitude and phase for a single sinusoid.  By
default the angular frequency is omega_0 but it can be any number or
symbol.

Phasors straddle the time and frequency domains and hence the
confusing phasor classes.

A phasor is described by an amplitude and phase.  There is also an
implicit frequency.  The amplitude and phase are usually functions of
frequency.  A phasor is useful to describe an AC voltage or current
signal.

A ratio of two phasors (of the same frequency) is no longer a
time-domain signal but a frequency domain quantity.  The frequency
dependence is explicit.  A phasor ratio is useful to describe an
immittance or transfer function.  These are frequency domain concepts.
A phasor ratio can be inferred from the Laplace domain by substituting
jomega (or jw) for s, where omega is the angular frequency of the phasor.

Lcapy considers V(jw) or I(jw) a generic phasor but V(3j) or I(3j) a
specific phasor.  Similarly, X(jw) is a generic phasor ratio and X(3j)
is a specific phasor ratio, where X denotes an immittance or transfer
function.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .expressionclasses import expressionclasses
from .acdc import ACChecker
from .domains import PhasorTimeDomain, PhasorFrequencyDomain
from .sym import j, omegasym
from .expr import expr
from .functions import sin, cos, exp, sqrt
from .expr import Expr
from .omegaexpr import omega
from sympy import Expr as symExpr
from warnings import warn


__all__ = ('phasor', )

# The phasor domain is different from the Fourier and Laplace domain
# since there is an implicit angular frequency.  This is only needed
# for voltage and current (vorrent) expressions.

# The phasor domain immittance can be found from the Laplace domain
# impedance/admittance by substituting s with j * omega.  In many
# cases the result is the same as the Fourier domain immittance but
# not always.  The most important case is the impedance of a
# capacitor.  In the Laplace domain it is 1/(s C); in the phasor
# domain it is 1/(j omega C); in the Fourier domain it is 1/(2 j omega
# C) + delta(omega) / (2 C).


class PhasorExpression(Expr):

    @property
    def omega(self):
        """Return angular frequency."""

        return self.assumptions.get('omega', None)

    @property
    def var(self):
        """Return underlying variable as sympy expression."""

        # TODO, think this through...
        var = self.omega
        if hasattr(var, 'expr'):
            raise ValueError('FIXME, omega is Lcapy symbol')

        # Handle things like 2 * pi * f
        if isinstance(var, symExpr) and var.is_Mul:
            return None

        return var

    def phasor(self, **assumptions):
        """Convert to phasor representation."""

        ass = self.assumptions.copy()
        assumptions = ass.merge(**assumptions)

        return self.__class__(self, **assumptions)

    def rms(self):
        """Return root mean square."""

        return abs(self) * sqrt(2) / 2

    def plot(self, wvector=None, **kwargs):
        """Plot polar diagram for a time-domain phasor or frequency response
        for a frequency-domain phasor.  For the latter, wvector
        specifies the angular frequencies.  If it is a tuple, it sets
        the angular frequency limits."""

        from .plot import plot_phasor, plot_angular_frequency

        if self.is_phasor_time_domain:
            return plot_phasor(self, **kwargs)

        if self.omega.is_constant():
            raise ValueError('Cannot plot at single frequency')

        return plot_angular_frequency(self, wvector, **kwargs)

    def bode_plot(self, fvector=None, unwrap=True, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Bode
        plot (but without the straight line approximations).  fvector
        specifies the frequencies.  If it is a tuple (m1, m2), it sets the
        frequency limits as (10**m1, 10**m2).

        `unwrap` controls phase unwrapping (default True).
        """

        from .plot import plot_bode

        if not self.is_phasor_frequency_domain:
            raise ValueError('Not frequency domain phasor: use plot()')

        result = self.fourier()
        return plot_bode(result, fvector, unwrap=unwrap, **kwargs)

    def _mul_compatible(self, x):

        if not hasattr(x, 'omega'):
            return True

        if self.omega == x.omega:
            return True

        raise ValueError('Incompatible phasor angular frequencies %s and %s' %
                         (self.omega, x.omega))

    def subs(self, *args, safe=False, **kwargs):
        """Substitute variables in expression, see sympy.subs for usage."""

        from .sym import fsym, pi

        # Handle omega -> 2 * pi * f conversions to Fourier domain.
        if len(args) == 1 and args[0] == 2 * pi * fsym:
            from .fexpr import FourierDomainExpression

            if self.var != omegasym:
                raise ValueError('Expecting omega for self.var')

            # Could check for expressions that are known to fail,
            # such as 1 / (j * omega).
            if not safe:
                warn("""
Converting to Fourier domain via phasor domain may not give correct answer.   It is safer to convert to time domain then to Fourier domain.""")

            return FourierDomainExpression(self.sympy.subs(omegasym, args[0].sympy)).as_quantity(self.quantity)

        return super(PhasorExpression, self).subs(*args, **kwargs)


class PhasorTimeDomainExpression(PhasorTimeDomain, PhasorExpression):
    """This is a phasor domain base class for voltages and currents."""

    is_phasor_domain = True
    is_phasor_time_domain = True

    def __init__(self, val, **assumptions):

        if isinstance(val, PhasorExpression):
            assumptions['omega'] = val.omega
        elif 'omega' not in assumptions or assumptions['omega'] is None:
            assumptions['omega'] = omegasym

        if hasattr(assumptions['omega'], 'expr'):
            assumptions['omega'] = assumptions['omega'].expr

        assumptions['ac'] = True
        super(PhasorExpression, self).__init__(val, **assumptions)

    def _class_by_quantity(self, quantity, domain=None):

        if quantity == 'undefined':
            return PhasorTimeDomainExpression

        return super(PhasorTimeDomainExpression, self)._class_by_quantity(quantity, domain)

    def as_expr(self):
        return PhasorTimeDomainExpression(self)

    @classmethod
    def from_time(cls, expr, omega=None, **assumptions):

        from .symbols import t

        if expr.is_admittance or expr.is_impedance or expr.is_transfer:
            warn('Should convert %s expression to Laplace-domain first.' %
                 expr.quantity)

        assumptions['ac'] = True

        if expr.is_transform_domain:
            warn('Converting %s-domain to time-domain first.' % expr.domain)
            expr = expr.time()

        check = ACChecker(expr, t)
        if not check.is_ac:
            raise ValueError(
                'Do not know how to convert %s to phasor.  Expecting an AC signal.' % expr)

        if omega is not None and check.omega != omega:
            raise ValueError('Expecting omega=%s, found omega=%s.' %
                             (omega, check.omega))

        if check.omega == 0:
            warn('DC phasor.')

        result = check.amp * exp(j * check.phase)
        assumptions['omega'] = check.omega
        assumptions['complex_signal'] = check.is_complex

        return cls.change(expr, result, domain='phasor', **assumptions)

    def time(self, **assumptions):
        """Convert to time domain representation."""
        from .symbols import t

        omega1 = self.omega
        if isinstance(omega1, Expr):
            # TODO: Fix inconsistency.  Sometimes omega is a symbol.
            omega1 = omega1.expr

        if self.is_complex_signal:
            result = self.expr * exp(j * omega1 * t)
        else:
            result = self.real.expr * \
                cos(omega1 * t) - self.imag.expr * sin(omega1 * t)

        return TimeDomainExpression(result).as_quantity(self.quantity)

    def laplace(self, **assumptions):
        """Convert to Laplace domain representation."""

        return self.time().laplace()

    def angular_fourier(self, **assumptions):
        """Angular Fourier transform."""

        if self.has(omega):
            warn('Expression contains omega, should substitute with a different symbol.')

        return self.time().angular_fourier()

    def fourier(self, **assumptions):
        """Fourier transform."""

        return self.time().fourier()


class PhasorFrequencyDomainExpression(PhasorFrequencyDomain, PhasorExpression):
    """This represents the ratio of two-phasors; for example
    an impedance, an admittance, or a transfer function."""

    is_phasor_frequency_domain = True
    is_transform_domain = True

    def __init__(self, val, **assumptions):

        if isinstance(val, PhasorExpression):
            assumptions['omega'] = val.omega
        elif 'omega' not in assumptions or assumptions['omega'] is None:
            assumptions['omega'] = omegasym

        if isinstance(val, Expr):
            ass = val.assumptions.copy()
            ass = ass.merge(**assumptions)
        else:
            ass = assumptions

        super(PhasorExpression, self).__init__(val, **ass)

    def _class_by_quantity(self, quantity, domain=None):

        if quantity == 'undefined':
            return PhasorFrequencyDomainExpression

        return super(PhasorFrequencyDomainExpression, self)._class_by_quantity(quantity, domain)

    @classmethod
    def from_laplace(cls, expr, **assumptions):

        from .sym import ssym

        ass = expr.assumptions.copy()
        ass = ass.merge(**assumptions)
        omega = ass.pop('omega', None)

        if omega is None:
            omega = omegasym
        elif hasattr(omega, 'expr'):
            omega = omega.expr

        if expr.is_voltage or expr.is_current:
            warn('Should convert %s expression to time-domain first.' %
                 expr.quantity)

        # Substitute jw for s
        result = expr.laplace(**ass)
        result2 = result.expr.replace(ssym, j * omegasym)

        quantity = expr.quantity
        if quantity == 'undefined':
            cls = PhasorFrequencyDomainExpression
        else:
            cls = expr._class_by_domain('phasor')

        ret = cls(result2, omega=omega, **ass)
        return ret

    def time(self, **assumptions):
        """Convert to time domain representation."""

        return self.laplace().time()

    def laplace(self, **assumptions):
        """Convert to Laplace domain representation."""

        from .sym import ssym
        from .sexpr import LaplaceDomainExpression

        ass = self.assumptions.copy()
        assumptions = ass.merge(**assumptions)

        omega = self.omega
        result = self.expr.replace(omega, ssym / j)
        return LaplaceDomainExpression(result, **assumptions).as_quantity(self.quantity)

    def angular_fourier(self, **assumptions):
        """Angular Fourier transform."""

        if self.var is omega:
            # TODO, add warning.
            pass

        return self.time().angular_fourier()

    def fourier(self, **assumptions):
        """Fourier transform."""

        return self.time().fourier()

    def as_expr(self):
        return PhasorFrequencyDomainExpression(self)


def phasor(arg, omega=None, **assumptions):
    """Create phasor.

    If arg has the form A * cos(w * t + phi) the phasor
    A * exp(j * phi) of angular frequency w is returned.

    If arg has the form A * sin(w * t + phi) the phasor
    A * exp(j * (phi - pi / 2)) of angular frequency w is returned.

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
        # Is this sensible?  It is probably better to have
        # a phasorratio function.  TODO add deprecation warning.
        return PhasorFrequencyDomainExpression(arg, omega=omega, **assumptions)


expressionclasses.register('phasor', PhasorTimeDomainExpression,
                           PhasorFrequencyDomainExpression,
                           ('voltage', 'current',
                            'voltagesquared', 'currentsquared', 'undefined'))

from .texpr import TimeDomainExpression  # nopep8
from .expr import Expr  # nopep8
