"""This module provides the FourierNoiseDomainExpression class to represent
f-domain (Fourier domain) noise expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""
from __future__ import division
from .sym import symsimplify
from .functions import sqrt
from .sym import pi, omegasym, fsym
from .state import state
from .domains import FourierNoiseDomain
from .expr import expr
from .noiseexpr import NoiseExpression
from .omegaexpr import omega, AngularFourierDomainExpression
from .fexpr import FourierDomainExpression
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
import sympy as sym
import numpy as np

class FourierNoiseDomainExpression(FourierNoiseDomain, NoiseExpression):
    """Frequency domain (one-sided) noise spectrum expression (amplitude
    spectral density).

    This characterises a wide-sense stationary, zero-mean Gaussian
    noise random process.

    When performing arithmetic on two FourierNoiseDomainExpression
    expressions it is assumed that they are uncorrelated unless they
    have the same nid (noise indentifier).  If the nid is not
    specified, a new one is created.

    Uncorrelated noise expressions are added in quadrature (on a power
    basis).  Thus (FourierNoiseDomainExpression(3) +
    FourierNoiseDomainExpression(4)).expr = 5 since 5 = sqrt(3**2 +
    4**2)

    FourierNoiseDomainExpression(3) != FourierNoiseDomainExpression(3)
    since they are different noise realisations albeit with the same
    properties.  However, FourierNoiseDomainExpression(3).expr ==
    FourierNoiseDomainExpression(3).expr.  Similarly,
    FourierNoiseDomainExpression(3, nid='n1') ==
    FourierNoiseDomainExpression(3, nid='n1') since they have the same
    noise identifier and thus have the same realisation.

    Caution: The sum of two noise expressions generates a noise
    expression with a new nid.  This can lead to unexpected results
    since noise expressions with different nids are assumed to be
    uncorrelated.  For example, consider:
    a = FourierNoiseDomainExpression(3);
    b = FourierNoiseDomainExpression(4)
    a + b - b gives sqrt(41) but a + b - a gives sqrt(34).

    This case is correctly handled by the SuperpositionVoltage and
    SuperpositionCurrent classes since each noise component is stored
    and considered separately.

    (SuperpositionVoltage(a) + SuperpositionVoltage(b) -
    SuperpositionVoltage(b)).n gives 3 as expected.

    """
    var = fsym

    def plot(self, fvector=None, **kwargs):
        """Plot frequency response at values specified by fvector.

        There are many plotting options, see matplotlib.pyplot.plot.

        For example:
            V.plot(fvector, log_frequency=True)
            V.real.plot(fvector, color='black')
            V.phase.plot(fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.
        """
        from .plot import plot_frequency

        return plot_frequency(self, fvector, **kwargs)

    def transform(self, arg, **assumptions):
        """Transform into a different domain."""

        arg = expr(arg)        
        if isinstance(arg, AngularFourierDomainExpression):
            result = self.subs(omega / (2 * pi))
            cls = self._class_by_quantity(self.quantity, 'angular fourier noise')
            return cls(result, nid=self.nid, **assumptions)
        elif isinstance(arg, FourierDomainExpression):
            result = self.subs(arg, **assumptions)
            cls = self._class_by_quantity(self.quantity, 'fourier noise')
            return cls(result, nid=self.nid, **assumptions)            

        return super(FourierNoiseDomainExpression, self).transform(arg, **assumptions)
    
    
class FourierNoiseDomainVoltage(VoltageMixin, FourierNoiseDomainExpression):
    """Voltage noise amplitude spectral density (units V/rtHz).
    This can be a function of linear frequency, f.  For example,
    to model an opamp voltage noise:

    v = FourierNoiseDomainVoltage(1e-8 / sqrt(f) + 8e-9)
    
    """

    quantity_label = 'Voltage noise spectral density'
    units = 'V/rtHz'


class FourierNoiseDomainCurrent(CurrentMixin, FourierNoiseDomainExpression):
    """Current noise amplitude spectral density (units A/rtHz).

    This can be a function of linear frequency, f.  For example,
    to model an opamp current noise:

    i = FourierNoiseDomainCurrent(3e-12 / sqrt(f) + 200e-15)
    """

    quantity_label = 'Current noise spectral density'
    units = 'A/rtHz'


from .expressionclasses import expressionclasses

expressionclasses.register('fourier noise', FourierNoiseDomainExpression, None,
                           ('voltage', 'current'))


from .fexpr import f
from .noiseomegaexpr import AngularFourierNoiseDomainExpression
