"""This module provides the FourierDomainExpression class to represent
f-domain (Fourier domain) expressions.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import FourierDomain
from .fourier import inverse_fourier_transform
from .expr import Expr, expr, expr_make
from .sym import fsym, ssym, tsym, pi
from .units import u as uu
from sympy import Integral, Expr as symExpr

class FourierDomainExpression(FourierDomain, Expr):

    """Fourier domain expression or symbol."""

    var = fsym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)        
        assumptions['real'] = True
        super(FourierDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr        
        if check and expr.find(ssym) != set() and not expr.has(Integral):
            raise ValueError(
                'f-domain expression %s cannot depend on s' % expr)
        if check and expr.find(tsym) != set() and not expr.has(Integral):
            raise ValueError(
                'f-domain expression %s cannot depend on t' % expr)

    def as_expr(self):
        return FourierDomainExpression(self)

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        from .symbols import omega
        
        result = self.subs(omega / (2 * pi))
        return result

    def inverse_fourier(self, evaluate=True, **assumptions):
        """Attempt inverse Fourier transform."""

        result = inverse_fourier_transform(self.expr, self.var, tsym, evaluate=evaluate)

        return self.change(result, 'time', units_scale=uu.Hz, **assumptions)

    def IFT(self, **assumptions):
        """Convert to t-domain.   This is an alias for inverse_fourier."""

        return self.inverse_fourier(**assumptions)    
    
    def time(self, **assumptions):
        return self.inverse_fourier(**assumptions)

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        from .symbols import omega
        
        result = self.subs(omega / (2 * pi))
        return result
    
    def laplace(self, **assumptions):
        """Determine one-side Laplace transform with 0- as the lower limit."""

        result = self.time(**assumptions).laplace()
        return result
    
    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        return self.time(**assumptions).phasor(**assumptions)        

    def plot(self, fvector=None, **kwargs):
        """Plot frequency response at values specified by fvector.  If fvector
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
        
        The plot axes are returned.  This is a tuple for magnitude/phase or
        real/imaginary plots.

        There are many plotting options, see lcapy.plot and
        matplotlib.pyplot.plot.

        For example:
            V.plot(fvector, log_frequency=True)
            V.real.plot(fvector, color='black')
            V.phase.plot(fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.

        """

        from .plot import plot_frequency
        return plot_frequency(self, fvector, **kwargs)


def fexpr(arg, **assumptions):
    """Create FourierDomainExpression object.  If `arg` is fsym return f"""

    if arg is fsym:
        return f
    return expr_make('fourier', arg, **assumptions)

from .expressionclasses import expressionclasses

classes = expressionclasses.register('fourier', FourierDomainExpression)

f = FourierDomainExpression('f')
f.units = uu.Hz
