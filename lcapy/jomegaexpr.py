"""This module provides the AngularFrequencyResponseDomainExpression class to
represent jomega-domain (angular frequency frequency domain) expressions.

Copyright 2022 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import AngularFrequencyResponseDomain
from .inverse_fourier import inverse_fourier_transform
from .expr import Expr, expr, expr_make
from .state import state, validate
from .sym import fsym, ssym, tsym, omegasym, j, pi
from .units import u as uu
from sympy import Expr as symExpr

__all__ = ()


class AngularFrequencyResponseDomainExpression(AngularFrequencyResponseDomain, Expr):
    """Frequency domain expression or symbol (angular frequency)."""

    var = omegasym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)

        super(AngularFrequencyResponseDomainExpression,
              self).__init__(val, **assumptions)

        expr = self.expr

        if check and not expr.has(omegasym):
            if expr.has(tsym):
                validate(state.t_in_jw,
                         'jomega-domain expression %s depends on t' % expr)
            if expr.has(ssym):
                validate(state.s_in_jw,
                         'jomega-domain expression % s depends on s' % expr)
            if expr.has(omegasym):
                validate(state.f_in_jw,
                         'jomega-domain expression %s depends on f' % expr)

    def as_expr(self):
        return AngularFrequencyResponseDomainExpression(self)

    def inverse_fourier(self, **assumptions):
        """Attempt inverse Fourier transform."""

        expr = self.subs(2 * pi * fsym)
        result = inverse_fourier_transform(expr.sympy, fsym, tsym)

        return self.change(result, 'time', units_scale=uu.Hz, **assumptions)

    def time(self, **assumptions):
        """Alias for inverse_fourier."""

        return self.inverse_fourier(**assumptions)

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""

        return self.laplace(**assumptions).angular_fourier()

    def fourier(self, **assumptions):
        """Convert to Fourier domain."""

        return self.laplace(**assumptions).fourier()

    def norm_fourier(self, **assumptions):
        """Convert to normalized Fourier domain."""

        return self.laplace(**assumptions).norm_fourier()

    def norm_angular_fourier(self, **assumptions):
        """Convert to normalized angular Fourier domain."""

        return self.laplace(**assumptions).norm_angular_fourier()

    def laplace(self, **assumptions):
        """Convert to Laplace domain."""

        from .symbols import s, j

        return self.subs(s / j)

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

    def bode_plot(self, wvector=None, unwrap=True, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Bode
        plot (but without the straight line approximations).  wvector
        specifies the angular frequencies.  If it is a tuple (f1, f2), it sets
        the frequency limits.   Since a logarithmic frequency scale is used,
        f1 must be greater than 0.

        `unwrap` controls phase unwrapping (default True).

        For more info, see `plot`.
        """

        from .plot import plot_angular_bode
        return plot_angular_bode(self, wvector, unwrap=unwrap, **kwargs)

    def nyquist_plot(self, wvector=None, log_frequency=True, **kwargs):
        """Plot frequency response as a Nyquist plot (imaginary part versus
        real part).  wvector specifies the angular frequencies.  If it
        is a tuple (f1, f2), it sets the frequency limits as (f1, f2).

        `npoints` set the number of plotted points.

        The unit circle is shown by default.  This can be disabled
        with `unitcircle=False`.

        """

        from .plot import plot_nyquist
        return plot_nyquist(self, wvector, log_frequency=log_frequency, **kwargs)

    def nichols_plot(self, wvector=None, log_frequency=True, **kwargs):
        """Plot frequency response as a Nichols plot (dB versus phase).
        wvector specifies the angular frequencies.  If it is a tuple
        (f1, f2), it sets the frequency limits as (f1, f2).

        `npoints` set the number of plotted points.

        """

        from .plot import plot_nichols
        return plot_nichols(self, wvector, log_frequency=log_frequency, **kwargs)


from .expressionclasses import expressionclasses  # nopep8

classes = expressionclasses.register('angular frequency response',
                                     AngularFrequencyResponseDomainExpression)

jomega = AngularFrequencyResponseDomainExpression('j * omega')
jomega.units = uu.rad / uu.s
