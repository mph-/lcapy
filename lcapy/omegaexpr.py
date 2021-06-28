"""This module provides the AngularFourierDomainExpression class to
represent omega-domain (angular frequency Fourier domain) expressions.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import AngularFourierDomain
from .inverse_fourier import inverse_fourier_transform
from .expr import Expr, expr, expr_make
from .sym import fsym, ssym, tsym, omegasym, omega0sym, j, pi
from .units import u as uu
from sympy import Expr as symExpr

__all__ = ('omegaexpr', )

class AngularFourierDomainExpression(AngularFourierDomain, Expr):
    """Fourier domain expression or symbol (angular frequency)."""

    var = omegasym

    def __init__(self, val, **assumptions):

        super(AngularFourierDomainExpression, self).__init__(val, **assumptions)

        if self.expr.has(ssym):
            raise ValueError(
                'omega-domain expression %s cannot depend on s' % self.expr)
        if self.expr.has(tsym):
            raise ValueError(
                'omega-domain expression %s cannot depend on t' % self.expr)

    def as_expr(self):
        return AngularFourierDomainExpression(self)
    
    def inverse_fourier(self, **assumptions):
        """Attempt inverse Fourier transform."""

        expr = self.subs(2 * pi * fsym)
        result = inverse_fourier_transform(expr, fsym, tsym)

        return self.change(result, 'time', units_scale=uu.Hz, **assumptions)

    def time(self, **assumptions):
        """Alias for inverse_fourier."""

        return self.inverse_fourier(**assumptions)

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""

        return self
    
    def fourier(self, **assumptions):
        """Convert to Fourier domain."""
        
        from .symbols import f
        
        result = self.subs(2 * pi * f)
        return result

    def norm_fourier(self, **assumptions):
        """Convert to normalized Fourier domain."""
        from .symbols import F
        from .dsym import dt
        
        result = self.subs(2 * pi * F / dt)
        return result

    def norm_angular_fourier(self, **assumptions):
        """Convert to normalized angular Fourier domain."""
        from .symbols import Omega
        from .dsym import dt
        
        result = self.subs(Omega / dt)
        return result        

    def laplace(self, **assumptions):
        """Convert to Laplace domain."""

        result = self.time(**assumptions).laplace()
        return result

    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        return self.time(**assumptions).phasor(**assumptions)

    def plot(self, wvector=None, **kwargs):
        """Plot angular frequency response at values specified by wvector.

        There are many plotting options, see matplotlib.pyplot.plot.

        For example:
            V.plot(fvector, log_frequency=True)
            V.real.plot(fvector, color='black')
            V.phase.plot(fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.

        `kwargs` include:
        `axes` - the plot axes to use otherwise a new figure is created
        `xlabel` - the x-axis label
        `ylabel` - the y-axis label
        `ylabel2` - the second y-axis label if needed, say for mag and phase
        `xscale` - the x-axis scaling, say for plotting as ms
        `yscale` - the y-axis scaling, say for plotting mV
        `norm` - use normalized frequency
        `dbmin` - the smallest value to plot in dB (default -120)
        in addition to those supported by the matplotlib plot command.

        The plot axes are returned.  This is a tuple for magnitude/phase or
        real/imaginary plots.
        """

        from .plot import plot_angular_frequency
        return plot_angular_frequency(self, wvector, **kwargs)

    def bode_plot(self, wvector=None, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Bode
        plot (but without the straight line approximations).  wvector
        specifies the angular frequencies.  If it is a tuple (f1, f2), it sets
        the frequency limits.   Since a logarithmic frequency scale is used,
        f1 must be greater than 0.

        For more info, see `plot`.
        """
        
        from .plot import plot_angular_bode
        return plot_angular_bode(self, wvector, **kwargs)

    def nyquist_plot(self, wvector=None, log_frequency=True, **kwargs):
        """Plot frequency response as a Nyquist plot (real part versus
        imaginary part).  wvector specifies the angular frequencies.  If it is
        a tuple (f1, f2), it sets the frequency limits as (f1, f2).

        `npoints` set the number of plotted points.

        The unit circle is shown by default.  This can be disabled with `unitcircle=False`.
        """        

        from .plot import plot_nyquist
        return plot_nyquist(self, wvector, log_frequency=log_frequency, **kwargs)                
    
        
def omegaexpr(arg, **assumptions):
    """Create AngularFourierDomainExpression object.  If `arg` is omegasym return omega"""

    if arg is omegasym:
        return omega
    if arg is omega0sym:
        return omega0
    return expr_make('angular fourier', arg, **assumptions)


from .expressionclasses import expressionclasses

classes = expressionclasses.register('angular fourier', AngularFourierDomainExpression)
        
omega = AngularFourierDomainExpression('omega')
omega.units = uu.rad / uu.s

# This represents a specific angular frequency and is assumed to be positive
omega0 = AngularFourierDomainExpression('omega_0', positive=True)
omega0.units = uu.rad / uu.s
