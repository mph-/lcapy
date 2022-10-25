"""This module provides the PhasorDomain and PhasorRatioDomain classes
 to represent phasors and ratios of phasors for AC analysis.

A phasor represents the amplitude and phase for a single sinusoid.  By
default the angular frequency is omega_0 but it can be any number or
symbol.

Phasors straddle the time and frequency domains and hence the
confusing phasor classes.

A phasor is described by an amplitude and phase.  There is also an
implicit frequency.  The amplitude and phase are usually functions of
frequency.  A phasor is useful to describe an AC voltage or current
signal.

Lcapy considers V(jw) or I(jw) a generic phasor but V(3j) or I(3j) a
specific phasor.  Similarly, X(jw) is a generic phasor ratio and X(3j)
is a specific phasor ratio, where X denotes an immittance or transfer
function.

A ratio of two phasors (of the same frequency) is not a phasor.  The
ratio of two generic phasors is a frequency domain quantity.  The
frequency dependence is explicit.  A phasor ratio is useful to
describe an immittance or transfer function.  These are frequency
domain concepts.  A phasor ratio can be inferred from the Laplace
domain by substituting jomega (or jw) for s, where omega is the
angular frequency of the phasor.

Copyright 2014--2022 Michael Hayes, UCECE

"""

from __future__ import division
from .expressionclasses import expressionclasses
from .acdc import ACChecker
from .domains import PhasorDomain, PhasorRatioDomain
from .sym import j, omegasym, fsym, pi
from .expr import expr
from .functions import sin, cos, exp, sqrt
from .expr import Expr
from .omegaexpr import omega
from .units import u as uu
from sympy import Expr as symExpr
from warnings import warn


__all__ = ('phasor', 'phasor_ratio')

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

        if self.omega == omegasym:
            return omegasym

        if self.omega == 2 * pi * fsym:
            return fsym

        return None

    def phasor(self, **assumptions):
        """Convert to phasor representation."""

        ass = self.assumptions.copy()
        assumptions = ass.merge(**assumptions)

        return self.__class__(self, **assumptions)

    def rms(self):
        """Return root mean square."""

        return abs(self) * sqrt(2) / 2

    def _compatible_phasors(self, x):

        if self.domain != x.domain:
            return False

        if not hasattr(x, 'omega'):
            return True

        if self.omega == x.omega:
            return True

        raise ValueError('Incompatible phasor angular frequencies %s and %s' %
                         (self.omega, x.omega))

    def _add_compatible_domains(self, x):

        return self._compatible_phasors(x)

    def _div_compatible_domains(self, x):

        if (self.is_constant_domain or x.is_constant_domain):
            return True

        return self._compatible_phasors(x)

    def _mul_compatible_domains(self, x):

        if (self.is_constant_domain or x.is_constant_domain):
            return True

        if self._compatible_phasors(x):
            return True

        if (self.is_phasor_ratio_domain and x.is_phasor_domain and
                self.omega == x.omega):
            return True
        if (self.is_phasor_domain and x.is_phasor_ratio_domain and
                self.omega == x.omega):
            return True

        # Allow phasor(x) * omega, etc.
        if (self.is_phasor_domain and x.is_angular_fourier_domain and
                self.omega == x.var):
            return True
        if (self.is_phasor_domain and x.is_fourier_domain and
                self.omega == 2 * pi * x.var):
            return True

        if (self.is_phasor_ratio_domain and x.is_angular_fourier_domain and
                self.omega == x.var):
            return True
        if (self.is_phasor_ratio_domain and x.is_fourier_domain and
                self.omega == 2 * pi * x.var):
            return True

        return False

    def subs(self, *args, safe=False, **kwargs):
        """Substitute variables in expression, see sympy.subs for usage."""

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


class PhasorDomainExpression(PhasorDomain, PhasorExpression):
    """This is a phasor domain base class for voltages and currents."""

    def __init__(self, val, **assumptions):

        if isinstance(val, PhasorExpression):
            assumptions['omega'] = val.omega
        elif 'omega' not in assumptions or assumptions['omega'] is None:
            assumptions['omega'] = omegasym

        if hasattr(assumptions['omega'], 'expr'):
            assumptions['omega'] = assumptions['omega'].expr

        assumptions['ac'] = True
        super(PhasorExpression, self).__init__(val, **assumptions)

    def _div_domain(self, x):

        if x.is_phasor_domain:
            return 'phasor ratio'

        return self.domain

    def _class_by_quantity(self, quantity, domain=None):

        if quantity == 'undefined':
            if domain == 'phasor ratio':
                return PhasorRatioDomainExpression
            else:
                return PhasorDomainExpression

        return super(PhasorDomainExpression, self)._class_by_quantity(quantity, domain)

    def as_expr(self):
        return PhasorDomainExpression(self)

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

    @classmethod
    def from_constant(cls, expr, omega=None, **assumptions):

        return cls(expr, omega=omega)

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

    def plot(self, **kwargs):
        """Plot polar diagram."""

        from .plot import plot_phasor

        return plot_phasor(self, **kwargs)


class PhasorRatioDomainExpression(PhasorRatioDomain, PhasorExpression):
    """This represents the ratio of two-phasors; for example
    an impedance, an admittance, or a transfer function."""

    is_phasor_ratio_domain = True
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

    def _mul_domain(self, x):

        if x.is_phasor_domain:
            return 'phasor'

        return self.domain

    def _class_by_quantity(self, quantity, domain=None):

        if quantity == 'undefined':
            return PhasorRatioDomainExpression

        return super(PhasorRatioDomainExpression, self)._class_by_quantity(quantity, domain)

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
        result2 = result.expr.replace(ssym, j * omega)

        quantity = expr.quantity
        if quantity == 'undefined':
            cls = PhasorRatioDomainExpression
        else:
            cls = expr._class_by_domain('phasor ratio')

        ret = cls(result2, omega=omega, **ass)
        return ret

    @property
    def domain_label(self):

        from .sym import fsym

        # TODO, remove this say by having PhasorRatioDomainExpression
        # and PhasorAngularFrequencyResponseDomainExpression classes.

        if self.var == fsym:
            return 'Frequency'
        else:
            return 'Angular frequency'

    @property
    def domain_units(self):

        from .sym import fsym
        from .units import u as uu

        # TODO, remove this say by having PhasorRatioDomainExpression
        # and PhasorAngularFrequencyResponseDomainExpression classes.

        if self.var == fsym:
            return uu.Hz
        else:
            return uu.rad / uu.s

    def time(self, **assumptions):
        """Convert to time domain representation."""

        return self.laplace().time()

    def laplace(self, **assumptions):
        """Convert to Laplace domain representation."""

        from .sym import ssym
        from .sexpr import LaplaceDomainExpression

        ass = self.assumptions.copy()
        assumptions = ass.merge(**assumptions)

        if self.var == omegasym:
            result = self.expr.replace(omegasym, ssym / j)
        elif self.var == fsym:
            result = self.expr.replace(fsym, ssym / (j * 2 * pi))
        else:
            raise ValueError(
                'Cannot convert PhasorRatioDomain to LaplaceDomain for variable %s' % self.var)

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
        return PhasorRatioDomainExpression(self)

    def plot(self, fvector=None, **kwargs):
        """Plot frequency response.

        If the dependent variable is omega, fvector specifies the angular
        frequencies.  If it is a tuple, it sets the angular frequency
        limits.

        If the dependent variable is f, fvector specifies the linear
        frequencies.  If it is a tuple, it sets the linear frequency
        limits."""

        from .plot import plot_angular_frequency, plot_frequency
        from .sym import fsym

        if self.omega.is_constant():
            raise ValueError('Cannot plot at single frequency')

        if self.var == fsym:
            return plot_frequency(self, fvector, **kwargs)

        return plot_angular_frequency(self, fvector, **kwargs)

    def bode_plot(self, fvector=None, unwrap=True, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Bode
        plot (but without the straight line approximations).  fvector
        specifies the frequencies.  If it is a tuple (m1, m2), it sets the
        frequency limits as (10**m1, 10**m2).

        `unwrap` controls phase unwrapping (default True).
        """

        from .plot import plot_bode

        return plot_bode(self, fvector, unwrap=unwrap, **kwargs)


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
        return PhasorDomainExpression.from_time(arg, omega=omega, **assumptions)
    else:
        # Expecting phasor (complex amplitude)
        return PhasorDomainExpression(arg, omega=omega, **assumptions)


def phasor_ratio(arg, omega=None, **assumptions):
    """Create phasor ratio.

    If arg is a constant C, the phasor ratio C is created with omega as the
    angular frequency (this defaults to 'omega' if not specified).

    """

    arg = expr(arg)

    return PhasorRatioDomainExpression(arg, omega=omega, **assumptions)


expressionclasses.register('phasor', PhasorDomainExpression,
                           None,
                           ('voltage', 'current',
                            'voltagesquared', 'currentsquared', 'power',
                            'undefined'))

# Handle voltage, current frequency response as phasor ratios
# so that voltage(1 / (s + 2)).subs(jw) works.
# Probably should add a frequency response domain.
expressionclasses.register('phasor ratio', PhasorRatioDomainExpression)

from .texpr import TimeDomainExpression  # nopep8
from .expr import Expr  # nopep8
