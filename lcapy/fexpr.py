"""This module provides the FourierDomainExpression class to represent
f-domain (Fourier domain) expressions.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import FourierDomain
from .inverse_fourier import inverse_fourier_transform
from .inverse_dtft import IDTFT
from .expr import Expr, expr, expr_make
from .state import state, validate
from .sym import fsym, ssym, tsym, omegasym, pi
from .sym import nsym, dt
from .units import u as uu
from .utils import factor_const
from sympy import Integral, Expr as symExpr


class FourierDomainExpression(FourierDomain, Expr):

    """Fourier domain expression or symbol."""

    var = fsym

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)
        assumptions['real'] = True
        super(FourierDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr

        if check and not expr.has(fsym):
            if expr.has(ssym):
                validate(state.s_in_f,
                         'f-domain expression % s depends on s' % expr)
            if expr.has(tsym):
                validate(state.t_in_f,
                         'f-domain expression %s depends on t' % expr)
            if expr.has(omegasym):
                validate(state.w_in_f,
                         'f-domain expression %s depends on omega' % expr)

    def as_expr(self):
        return FourierDomainExpression(self)

    def inverse_fourier(self, evaluate=True, **assumptions):
        """Attempt inverse Fourier transform."""

        result = inverse_fourier_transform(
            self.expr, self.var, tsym, evaluate=evaluate)

        return self.change(result, 'time', units_scale=uu.Hz, **assumptions)

    def IFT(self, evaluate=True, **assumptions):
        """Convert to time domain.  This is an alias for inverse_fourier."""

        return self.inverse_fourier(evaluate=evaluate, **assumptions)

    def IDTFT(self, evaluate=True, **assumptions):
        """Convert to discrete-time domain using inverse discrete-time
        Fourier transform."""

        result = IDTFT(self.expr, self.var, nsym, evaluate=evaluate)

        return self.change(result, 'discrete time', units_scale=uu.Hz,
                           **assumptions)

    def time(self, **assumptions):
        return self.inverse_fourier(**assumptions)

    def norm_fourier(self, **assumptions):
        """Convert to normalized Fourier domain."""
        from .symbols import F
        from .sym import dt

        result = self.subs(F / dt)
        return result

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        from .symbols import omega

        result = self.subs(omega / (2 * pi))
        return result

    def norm_angular_fourier(self, **assumptions):
        """Convert to normalized angular Fourier domain."""
        from .symbols import Omega
        from .sym import dt

        result = self.subs(Omega / (2 * pi * dt))
        return result

    def frequency_response(self, **assumptions):
        """Convert to frequency response domain."""

        return self.laplace(**assumptions).frequency_response()

    def angular_frequency_response(self, **assumptions):
        """Convert to angular frequency response domain."""

        return self.laplace(**assumptions).angular_frequency_response()

    def laplace(self, **assumptions):
        """Determine one-sided Laplace transform with 0- as the lower limit."""

        result = self.time(**assumptions).laplace()
        return result

    def phasor(self, **assumptions):
        """Convert to phasor domain."""

        return self.time(**assumptions).phasor(**assumptions)

    def plot(self, fvector=None, plot_type=None, **kwargs):
        """Plot frequency response at values specified by `fvector`.

        If `fvector` is a tuple, this sets the frequency limits.

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
        `norm` - use normalized frequency
        `dbmin` - the smallest value to plot in dB (default -120)
        `plot_deltas` - plot Dirac deltas as arrows
        in addition to those supported by the matplotlib plot command.

        The plot axes are returned.  This is a tuple for magnitude/phase or
        real/imaginary plots.

        There are many plotting options, see lcapy.plot and
        matplotlib.pyplot.plot.

        For example:
            V.plot(fvector, log_frequency=True)
            V.real.plot(fvector, color='black')
            V.phase.plot(fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of
        magnitude (dB) and phase.

        """

        from .plot import plot_frequency
        return plot_frequency(self, fvector, plot_type=plot_type, **kwargs)

    def bode_plot(self, fvector=None, unwrap=True, **kwargs):
        """Plot frequency response for a frequency-domain phasor as a Bode
        plot (but without the straight line approximations).  fvector
        specifies the frequencies.  If it is a tuple (f1, f2), it sets
        the frequency limits.   Since a logarithmic frequency scale is used,
        f1 must be greater than 0.

        `unwrap` controls phase unwrapping (default True).

        For more info, see `plot`.
        """

        from .plot import plot_bode
        return plot_bode(self, fvector, unwrap=unwrap, **kwargs)

    def nyquist_plot(self, fvector=None, log_frequency=True, **kwargs):
        """Plot frequency response as a Nyquist plot (imaginary part versus
        real part).  fvector specifies the frequencies.  If it is
        a tuple (f1, f2), it sets the frequency limits as (f1, f2).

        `npoints` set the number of plotted points.

        The unit circle is shown by default.  This can be disabled with `unitcircle=False`.
        """

        from .plot import plot_nyquist
        return plot_nyquist(self, fvector, log_frequency=log_frequency, **kwargs)

    def nichols_plot(self, fvector=None, log_frequency=True, **kwargs):
        """Plot frequency response as a Nichols plot (dB part versus
        phase).  fvector specifies the frequencies.  If it is
        a tuple (f1, f2), it sets the frequency limits as (f1, f2).

        `npoints` set the number of plotted points.
        """

        from .plot import plot_nichols
        return plot_nichols(self, fvector, log_frequency=log_frequency, **kwargs)


def fexpr(arg, **assumptions):
    """Create FourierDomainExpression object.  If `arg` is fsym return f"""

    if arg is fsym:
        return f
    return expr_make('fourier', arg, **assumptions)


from .expressionclasses import expressionclasses  # nopep8

classes = expressionclasses.register('fourier', FourierDomainExpression)

f = FourierDomainExpression('f')
f.units = uu.Hz
