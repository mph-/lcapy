"""This module provides the noisefExpr class to represent
f-domain (Fourier domain) noise expressions.

Copyright 2014--2019 Michael Hayes, UCECE

"""
from __future__ import division
from .sym import symsimplify
from .functions import sqrt
from .sym import pi, omegasym, fsym
from .state import state
from .expr import expr
from .noiseexpr import noiseExpr
from .omegaexpr import omega, omegaExpr
from .fexpr import fExpr
import sympy as sym
import numpy as np

class noisefExpr(noiseExpr):
    """Frequency domain (one-sided) noise spectrum expression (amplitude
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

    This case is correctly handled by the Super class since each noise
    component is stored and considered separately.

    (Voltage(a) + Voltage(b) - Voltage(b)).n gives 3 as expected.

    """
    one_sided = True
    var = fsym

    domain_name = 'Frequency'
    domain_units = 'Hz'    

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
        if isinstance(arg, omegaExpr):
            result = self.subs(omega / (2 * pi))
            return result.subs(arg, **assumptions)            
        elif isinstance(arg, fExpr):
            return self.subs(arg, **assumptions)

        return super(noisefExpr, self).transform(arg, **assumptions)
    

class Vfnoisy(noisefExpr):
    """Voltage noise amplitude spectral density (units V/rtHz).
    This can be a function of linear frequency, f.  For example,
    to model an opamp voltage noise:

    v = Vfnoisy(1e-8 / sqrt(f) + 8e-9)
    
    """

    quantity = 'Voltage noise spectral density'
    units = 'V/rtHz'

    
    def __init__(self, val, **assumptions):

        super(Vfnoisy, self).__init__(val, **assumptions)
        # FIXME
        self._fourier_conjugate_class = Vt
        self._subs_classes = {omegaExpr: Vnoisy, fExpr: Vfnoisy}


class Ifnoisy(noisefExpr):
    """Current noise amplitude spectral density (units A/rtHz).

    This can be a function of linear frequency, f.  For example,
    to model an opamp current noise:

    i = Ifnoisy(3e-12 / sqrt(f) + 200e-15)
    """

    quantity = 'Current noise spectral density'
    units = 'A/rtHz'

    def __init__(self, val, **assumptions):

        super(Ifnoisy, self).__init__(val, **assumptions)
        # FIXME
        self._fourier_conjugate_class = It
        self._subs_classes = {omegaExpr: Inoisy, fExpr: Ifnoisy}

    
from .texpr import It, Vt
from .fexpr import f
from .noiseomegaexpr import Vnoisy, Inoisy
