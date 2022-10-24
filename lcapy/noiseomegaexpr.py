"""This module provides the AngularFourierNoiseDomainExpression class to represent omega-domain (angular Fourier domain) noise expressions.

Copyright 2014--2022 Michael Hayes, UCECE
"""

from __future__ import division
from .sym import symsimplify
from .functions import sqrt
from .sym import pi, omegasym, fsym
from .state import state
from .domains import AngularFourierNoiseDomain
from .expr import expr
from .noiseexpr import NoiseExpression
from .fexpr import f, FourierDomainExpression
from .omegaexpr import AngularFourierDomainExpression
import sympy as sym


class AngularFourierNoiseDomainExpression(AngularFourierNoiseDomain, NoiseExpression):
    """Angular frequency domain (one-sided) noise spectrum expression (amplitude
    spectral density).

    This characterises a wide-sense stationary, zero-mean Gaussian
    noise random process.

    When performing arithmetic on two AngularFourierNoiseDomainExpression
    expressions it is assumed that they are uncorrelated unless they have
    the same nid (noise indentifier).  If the nid is not specified, a new one is created.

    Uncorrelated noise expressions are added in quadrature (on a power
    basis).  Thus (AngularFourierNoiseDomainExpression(3) +
    AngularFourierNoiseDomainExpression(4)).expr = 5 since 5 =
    sqrt(3**2 + 4**2)

    AngularFourierNoiseDomainExpression(3) !=
    AngularFourierNoiseDomainExpression(3) since they are different
    noise realisations albeit with the same properties.  However,
    AngularFourierNoiseDomainExpression(3).expr ==
    AngularFourierNoiseDomainExpression(3).expr.  Similarly,
    AngularFourierNoiseDomainExpression(3, nid='n1') ==
    AngularFourierNoiseDomainExpression(3, nid='n1') since they have
    the same noise identifier and thus have the same realisation.

    Caution: The sum of two noise expressions generates a noise
    expression with a new nid.  This can lead to unexpected results
    since noise expressions with different nids are assumed to be
    uncorrelated.  For example, consider:
    a = AngularFourierNoiseDomainExpression(3);
    b = AngularFourierNoiseDomainExpression(4)
    a + b - b gives sqrt(41) and a + b - a gives sqrt(34).

    This case is correctly handled by the SuperpositionVoltage and
    SuperpositionCurrent classes since each noise component is stored
    and considered separately.

    (SuperpositionVoltage(a) + SuperpositionVoltage(b) - SuperpositionVoltage(b)).n gives 3 as expected.

    """
    var = omegasym

    def plot(self, omegavector=None, **kwargs):
        """Plot frequency response at values specified by omegavector.

        There are many plotting options, see matplotlib.pyplot.plot.

        For example:
            V.plot(omegavector, log_frequency=True)
            V.real.plot(omegavector, color='black')
            V.phase.plot(omegavector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.
        """
        from .plot import plot_angular_frequency

        return plot_angular_frequency(self, omegavector, **kwargs)

    def transform(self, arg, **assumptions):
        """Transform into a different domain."""

        arg = expr(arg)
        if isinstance(arg, FourierDomainExpression):
            result = self.subs(f * 2 * pi)
            cls = self._class_by_quantity(self.quantity, 'fourier noise')
            return cls(result, nid=self.nid, **assumptions)
        elif isinstance(arg, AngularFourierDomainExpression):
            result = self.subs(arg, **assumptions)
            cls = self._class_by_quantity(
                self.quantity, 'angular fourier noise')
            return cls(result, nid=self.nid, **assumptions)

        return super(AngularFourierNoiseDomainExpression, self).transform(arg, **assumptions)


class AngularFourierNoiseDomainVoltage(AngularFourierNoiseDomainExpression):
    """Voltage noise amplitude spectral density (units V/rtrad/s).
    This can be a function of angular frequency, omega.  For example,
    to model an opamp voltage noise:

    v = AngularFourierNoiseDomainVoltage(1e-8 / sqrt(omega) + 8e-9)

    """

    quantity_label = 'Voltage noise spectral density'
    units = 'V/rtrad/s'


class AngularFourierNoiseDomainCurrent(AngularFourierNoiseDomainExpression):
    """Current noise amplitude spectral density (units A/rtrad/s).

    This can be a function of angular frequency, omega.  For example,
    to model an opamp current noise:

    i = AngularFourierNoiseDomainCurrent(3e-12 / sqrt(omega) + 200e-15)
    """

    quantity_label = 'Current noise spectral density'
    units = 'A/rtrad/s'


from .expressionclasses import expressionclasses  # nopep8

classes = expressionclasses.register('angular fourier noise',
                                     AngularFourierNoiseDomainExpression, None,
                                     ('voltage', 'current'))
AngularFourierNoiseDomainVoltage = classes['voltage']
AngularFourierNoiseDomainCurrent = classes['current']

from .omegaexpr import omega  # nopep8
from .noisefexpr import FourierNoiseDomainExpression  # nopep8
