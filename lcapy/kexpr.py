"""This module provides the DiscreteFourierDomainExpression class to
 represent k-domain (discrete Fourier domain) expressions.

Copyright 2020--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import DiscreteFourierDomain
from .fourier import inverse_fourier_transform
from .functions import exp
from .sym import j, oo, pi
from .seqexpr import SequenceExpression
from .dsym import nsym, ksym, zsym
from .dft import IDFT
from sympy import Sum

__all__ = ('kexpr', )

class DiscreteFourierDomainExpression(DiscreteFourierDomain, SequenceExpression):
    """Discrete Fourier domain expression or symbol."""

    var = ksym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)
        assumptions['real'] = True
        
        super(DiscreteFourierDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr
        if check and expr.find(zsym) != set():
            raise ValueError(
                'k-domain expression %s cannot depend on z' % expr)
        if check and expr.find(nsym) != set() and not expr.has(Sum):
            raise ValueError(
                'k-domain expression %s cannot depend on n' % expr)

    def as_expr(self):
        return DiscreteFourierDomainExpression(self)

    def plot(self, kvector=None, **kwargs):
        """Plot frequency response at values specified by kvector.  If kvector
        is a tuple, this sets the frequency limits.

        kwargs include:
        axes - the plot axes to use otherwise a new figure is created
        xlabel - the x-axis label
        ylabel - the y-axis label
        ylabel2 - the second y-axis label if needed, say for mag and phase
        xscale - the x-axis scaling, say for plotting as ms
        yscale - the y-axis scaling, say for plotting mV
        plot_type -  'dB_phase', 'mag-phase', 'real-imag', 'mag', 'phase',
        'real', or 'imag'
        in addition to those supported by the matplotlib plot command.
        
        The plot axes are returned.

        There are many plotting options, see lcapy.plot and
        matplotlib.pyplot.plot.

        For example:
            V.plot(kvector, log_frequency=True)
            V.real.plot(kvector, color='black')
            V.phase.plot(kvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.

        """

        from .plot import plot_frequency
        return plot_frequency(self, kvector, **kwargs)

    def IDFT(self, N=None, evaluate=True):

        from .nexpr import n

        if N is None:
            from .sym import sympify
            
            N = sympify('N')

        result = IDFT(self.expr, ksym, nsym, N, evaluate=evaluate)
        return self.change(result, domain='discrete time')
    
    def ZT(self, **assumptions):
        return self.IDFT().ZT(**assumptions)

    
def kexpr(arg, **assumptions):
    """Create kExpr object.  If `arg` is ksym return k"""

    if arg is ksym:
        return k
    return DiscreteFourierDomainExpression(arg, **assumptions)


from .expressionclasses import expressionclasses

classs = expressionclasses.register('discrete fourier', DiscreteFourierDomainExpression)

k = DiscreteFourierDomainExpression('k')
