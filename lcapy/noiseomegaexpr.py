"""This module provides the AngularFourierDomainNoiseExpression class to represent
omega-domain (angular Fourier domain) noise expressions.

Copyright 2014--2020 Michael Hayes, UCECE
"""

from __future__ import division
from .sym import symsimplify
from .functions import sqrt
from .sym import pi, omegasym, fsym
from .state import state
from .expr import expr
from .noiseexpr import NoiseExpression
from .fexpr import f, FourierDomainExpression
from .omegaexpr import AngularFourierDomainExpression
import sympy as sym
import numpy as np

class AngularFourierDomainNoiseExpression(NoiseExpression):
    """Angular frequency domain (one-sided) noise spectrum expression (amplitude
    spectral density).

    This characterises a zero-mean Gaussian noise process.

    When performing arithmetic on two noiseExpr expressions it is
    assumed that they are uncorrelated unless they have the same nid
    (noise indentifier).  If the nid is not specified, a new one is
    created.

    Uncorrelated noise expressions are added in quadrature (on a power
    basis).  Thus (Vnoisy(3) + Vnoisy(4)).expr = 5 since 5 = sqrt(3**2 + 4**2)

    Vnoisy(3) != Vnoisy(3) since they are different noise realisations albeit
    with the same properties.  However, Vnoisy(3).expr == Vnoisy(3).expr.
    Similarly, Vnoisy(3, nid='n1') == Vnoisy(3, nid='n1') since they have the
    same noise identifier and thus have the same realisation.

    Caution: The sum of two noise expressions generates a noise
    expression with a new nid.  This can lead to unexpected results
    since noise expressions with different nids are assumed to be
    uncorrelated.  For example, consider:
    a = Vnoisy(3); b = Vnoisy(4)
    a + b - b gives sqrt(41) and  a + b - a gives sqrt(34).

    This case is correctly handled by the SuperpositionVoltage and
    SuperpositionCurrent classes since each noise component is stored
    and considered separately.

    (SuperpositionVoltage(a) + SuperpositionVoltage(b) - SuperpositionVoltage(b)).n gives 3 as expected.

    """
    one_sided = True
    var = omegasym

    domain_label = 'Angular frequency'
    domain_units = 'rad/s'        


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
            return result.subs(arg, **assumptions)
        elif isinstance(arg, AngularFourierDomainExpression):
            return self.subs(arg, **assumptions)        

        return super(AngularFourierDomainNoiseExpression, self).transform(arg, **assumptions)
    

class AngularFourierDomainNoiseVoltage(AngularFourierDomainNoiseExpression):
    """Voltage noise amplitude spectral density (units V/rtrad/s).
    This can be a function of angular frequency, omega.  For example,
    to model an opamp voltage noise:

    v = AngularFourierDomainNoiseVoltage(1e-8 / sqrt(omega) + 8e-9)
    
    """

    quantity_label = 'Voltage noise spectral density'
    units = 'V/rtrad/s'

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainNoiseVoltage, self).__init__(val, **assumptions)
        self._subs_classes = {FourierDomainExpression: FourierDomainNoiseVoltage, AngularFourierDomainExpression: AngularFourierDomainNoiseVoltage}


class AngularFourierDomainNoiseCurrent(AngularFourierDomainNoiseExpression):
    """Current noise amplitude spectral density (units A/rtrad/s).

    This can be a function of angular frequency, omega.  For example,
    to model an opamp current noise:

    i = AngularFourierDomainNoiseCurrent(3e-12 / sqrt(omega) + 200e-15)
    """

    quantity_label = 'Current noise spectral density'
    units = 'A/rtrad/s'

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainNoiseCurrent, self).__init__(val, **assumptions)
        self._subs_classes = {FourierDomainExpression: FourierDomainNoiseCurrent, AngularFourierDomainExpression: AngularFourierDomainNoiseCurrent}
    

from .texpr import TimeDomainCurrent, TimeDomainVoltage
from .omegaexpr import omega
from .noisefexpr import FourierDomainNoiseVoltage, FourierDomainNoiseCurrent
