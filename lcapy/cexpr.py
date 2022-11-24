"""This module provides the ConstantDomainExpression class to
represent constant expressions.

Note there are two types:
1. ConstantTimeDomainExpression for signals
2. ConstantFrequencyResponseDomainExpression for transfer functions and immitances.

Note, impedance(3)(s) gives 3 but voltage(3)(s) gives 3 / s.
Similarly, voltage(3)(t) gives 3 but impedance(3)(t) gives 3 * delta)(t)

Copyright 2014--2022 Michael Hayes, UCECE

"""

from .expr import Expr, expr_make
from .domains import ConstantDomain
from .sym import symbols_find
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
from .admittancemixin import AdmittanceMixin
from .impedancemixin import ImpedanceMixin
from .transfermixin import TransferMixin


class ConstantDomainExpression(ConstantDomain, Expr):
    """Constant real expression or symbol.

    If symbols in the expression are known to be negative, use
    ConstantDomainExpression(expr, positive=False)

    """

    def __init__(self, val, **assumptions):

        symbols = symbols_find(val)
        for symbol in ('s', 'omega', 't', 'f'):
            if symbol in symbols:
                raise ValueError(
                    'Constant expression %s cannot depend on %s' % (val, symbol))

        assumptions['dc'] = True
        super(ConstantDomainExpression, self).__init__(val, **assumptions)

    def as_expr(self):
        return ConstantDomainExpression(self)

    def rms(self):
        """Return root mean square."""
        return self

    def fourier(self):
        """Convert to Fourier domain representation."""

        return self.time().fourier()

    def angular_fourier(self):
        """Convert to angular Fourier domain representation."""

        return self.time().angular_fourier()

    def frequency_response(self):
        """Convert to  frequency response domain representation."""

        return self.change(self, domain='angular frequency response')

    def angular_frequency_response(self):
        """Convert to angular frequency response domain representation."""

        return self.change(self, domain='angular frequency response')

    def canonical(self, factor_const=True):
        # Minor optimisation
        return self

    def phasor(self, omega=None, **assumptions):

        from .phasor import PhasorDomainExpression

        return PhasorDomainExpression.from_constant(self, omega=omega,
                                                    **assumptions)

    def phasor_ratio(self, omega=None, **assumptions):

        return self.as_laplace(**assumptions).phasor_ratio(omega=omega,
                                                           **assumptions)

    def time(self, **assumptions):
        """Convert to time domain representation."""

        return self.change(self, domain='time', **assumptions)

    def laplace(self, **assumptions):
        """Convert to Laplace domain representation."""

        return self.time().laplace(**assumptions)
        # return self.change(self, domain='laplace', **assumptions)

    def as_laplace(self, **assumptions):
        """Convert to Laplace domain representation."""

        return self.change(self, domain='laplace', **assumptions)

    def response(self, xvector, tvector):
        """Evaluate response to input signal `xvector` at times
        `tvector`.  This returns a NumPy array."""

        from .sexpr import s
        # This can be optimized!
        return self(s).response(xvector, tvector)


class ConstantTimeDomainExpression(ConstantDomainExpression):
    """This represents constant voltage, current, voltage-squared, and
    current-squared signals.

    Note, impedance(3)(s) gives 3 but voltage(3)(s) gives 3 / s.

    """

    def _class_by_quantity(self, quantity, domain=None):

        if quantity == 'undefined' and domain is None:
            return ConstantTimeDomainExpression

        return super(ConstantTimeDomainExpression, self)._class_by_quantity(quantity, domain)

    def phasor_ratio(self, omega=None, **assumptions):

        raise ValueError(
            'Cannot convert %s(%s) to phasor_ratio' %
            (self.__class__.__name__, self))


class ConstantFrequencyResponseDomainExpression(ConstantDomainExpression):
    """This represents constant impedance, admittance, and transfer
    functions.

    Note, impedance(3)(s) gives 3 but voltage(3)(s) gives 3 / s.
    """

    def _class_by_quantity(self, quantity, domain=None):

        if quantity == 'undefined' and domain is None:
            return ConstantFrequencyResponseDomainExpression

        return super(ConstantFrequencyResponseDomainExpression, self)._class_by_quantity(quantity, domain)

    def phasor(self, omega=None, **assumptions):

        raise ValueError(
            'Cannot convert %s(%s) to phasor' %
            (self.__class__.__name__, self))

    def laplace(self, **assumptions):
        return self.change(self, domain='laplace', **assumptions)

    def time(self, **assumptions):
        """Convert to time domain."""

        return self.laplace(**assumptions).time(**assumptions)


def cexpr(arg, frequency=False, **assumptions):
    """Create Lcapy constant expression from `arg`.

    By default, `arg` is assumed to be positive.  If symbols in the
    `arg` are known to be negative, use `cexpr(arg, positive=False)`.

    """
    try:
        quantity = arg.quantity
    except:
        quantity = 'undefined'

    if quantity == 'undefined':

        if frequency:
            return ConstantFrequencyResponseDomainExpression(arg, **assumptions)

        # Sit on the fence rather than choosing ConstantTimeDomainExpression
        return ConstantDomainExpression(arg, **assumptions)

    return expr_make('constant', arg, **assumptions)


from .expressionclasses import expressionclasses  # nopep8

expressionclasses.register('constant', ConstantTimeDomainExpression,
                           ConstantFrequencyResponseDomainExpression,
                           ('voltage', 'current', 'voltagesquared',
                            'currentsquared', 'undefined'))
