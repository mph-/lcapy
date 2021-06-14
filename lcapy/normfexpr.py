"""This module provides the NormFourierDomainExpression class
to represent F-domain (normalised  Fourier domain)
expressions.

Copyright 2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import NormFourierDomain
from .inverse_fourier import inverse_fourier_transform
from .inverse_dtft import IDTFT
from .expr import Expr, expr, expr_make
from .fexpr import f
from .sym import ssym, tsym, fsym, pi
from .dsym import Fsym, nsym, dt
from .units import u as uu
from .utils import factor_const, remove_images
from sympy import Integral, Expr as symExpr

class NormFourierDomainExpression(NormFourierDomain, Expr):

    """Normalised  Fourier domain expression or symbol."""

    var = Fsym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)        
        assumptions['real'] = True
        super(NormFourierDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr        
        if check and expr.find(ssym) != set() and not expr.has(Integral):

            raise ValueError(
                'F-domain expression %s cannot depend on s' % expr)
        if check and expr.find(tsym) != set() and not expr.has(Integral):
            raise ValueError(
                'F-domain expression %s cannot depend on t' % expr)

    def as_expr(self):
        return NormFourierDomainExpression(self)

    def inverse_fourier(self, evaluate=True, **assumptions):
        """Attempt inverse Fourier transform."""

        expr = self.subs(2 * pi * fsym * dt)        
        result = inverse_fourier_transform(expr, fsym, tsym, evaluate=evaluate)

        return self.change(result, 'time', units_scale=uu.Hz, **assumptions)

    def IFT(self, evaluate=True, **assumptions):
        """Convert to time domain.  This is an alias for inverse_fourier."""

        return self.inverse_fourier(evaluate=evaluate, **assumptions)

    def IDTFT(self, var=None, evaluate=True, **assumptions):
        """Convert to discrete-time domain using inverse discrete-time
        Fourier transform."""

        foo = self.subs(2 * pi * f * dt)
        result = IDTFT(foo.expr, self.var, nsym, evaluate=evaluate)

        return self.change(result, 'discrete time', units_scale=uu.Hz,
                           **assumptions)

    def time(self, **assumptions):
        return self.inverse_fourier(**assumptions)

    def fourier(self, **assumptions):
        """Convert to Fourier domain."""
        from .symbols import f
        
        result = self.subs(f * dt)
        return result

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        from .symbols import omega
        
        result = self.subs(omega / dt)
        return result

    def norm_fourier(self, **assumptions):
        """Convert to normalised Fourier domain."""
        from .symbols import F
        
        result = self.subs(F / dt)
        return result

    def norm_angular_fourier(self, **assumptions):
        """Convert to normalised angular Fourier domain."""

        return self
    
    def laplace(self, **assumptions):
        """Determine one-side Laplace transform with 0- as the lower limit."""

        result = self.time(**assumptions).laplace()
        return result
    
    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        return self.time(**assumptions).phasor(**assumptions)        

    def plot(self, Fvector=None, plot_type=None, **kwargs):
        """Plot frequency response at values specified by `Fvector`.

        If `Fvector` is a tuple, this sets the frequency limits.

        `plot_type` - 'dB-phase', 'dB-phase-degrees', 'mag-phase',
        'mag-phase-degrees', 'real-imag', 'mag', 'phase',
        'phase-degrees', 'real', or 'imag'.

        The default `plot_type` for complex data is `dB-phase`.

        `kwargs` include:
        `axes` - the plot axes to use otherwise a new figure is created
        `xlabel` - the x-axis label
        `ylabel` - the y-axis label
        `ylabel2` - the second y-axis label if needed, say for mag and phase
        `xscale` - the x-axis scaling, say for plotting as ms
        `yscale` - the y-axis scaling, say for plotting mV
        `norm` - use normalised frequency
        in addition to those supported by the matplotlib plot command.
        
        The plot axes are returned.  This is a tuple for magnitude/phase or
        real/imaginary plots.

        There are many plotting options, see lcapy.plot and
        matplotlib.pyplot.plot.

        For example:
            V.plot(Fvector, log_frequency=True)
            V.real.plot(Fvector, color='black')
            V.phase.plot(Fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of
        magnitude (dB) and phase.

        """

        from .plot import plot_frequency
        return plot_frequency(self, Fvector, plot_type=plot_type,
                              norm=True, **kwargs)

    def remove_images(self):
        """Remove all the spectral images resulting from a DTFT.
        
        For example,

        >>> x = Sum(DiracDelta(F - m * 2 * pi), (m, -oo, oo))
        >>> x.remove_images()
        DiracDelta(F)
        """

        result = remove_images(self.expr, self.var, dt)
        return self.__class__(result, **self.assumptions)

    
def Fexpr(arg, **assumptions):
    """Create NormFourierDomainExpression object. 
    If `arg` is Fsym return F"""

    if arg is Fsym:
        return F
    return expr_make('fourier', arg, **assumptions)

from .expressionclasses import expressionclasses

classes = expressionclasses.register('norm fourier',
                                     NormFourierDomainExpression)

F = NormFourierDomainExpression('F')
F.units = uu.rad / uu.rad

