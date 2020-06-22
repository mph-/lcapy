"""This module provides the fExpr class to represent f-domain (Fourier
domain) expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .fourier import inverse_fourier_transform
from .expr import Expr, expr
from .sym import fsym, ssym, tsym, pi
from sympy import Integral

class fExpr(Expr):

    """Fourier domain expression or symbol."""

    var = fsym
    domain_name = 'Frequency'
    domain_units = 'Hz'

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)        
        assumptions['real'] = True
        super(fExpr, self).__init__(val, **assumptions)
        # Define when class defined.
        self._fourier_conjugate_class = tExpr

        expr = self.expr        
        if check and expr.find(ssym) != set() and not expr.has(Integral):
            raise ValueError(
                'f-domain expression %s cannot depend on s' % expr)
        if check and expr.find(tsym) != set() and not expr.has(Integral):
            raise ValueError(
                'f-domain expression %s cannot depend on t' % expr)                            
    def inverse_fourier(self, evaluate=True, **assumptions):
        """Attempt inverse Fourier transform."""

        result = inverse_fourier_transform(self.expr, self.var, tsym, evaluate=evaluate)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        else:
            result = tExpr(result)
        return result

    def IFT(self, **assumptions):
        """Convert to t-domain.   This is an alias for inverse_fourier."""

        return self.inverse_fourier(**assumptions)    
    
    def time(self, **assumptions):
        return self.inverse_fourier()
    
    def laplace(self, **assumptions):
        """Determine one-side Laplace transform with 0- as the lower limit."""

        return self.time(**assumptions).laplace()
    
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
        
        The plot axes are returned.

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

    def transform(self, arg, **assumptions):
        """Transform into a different domain."""

        from .omegaexpr import omegaExpr, omega

        arg = expr(arg)        
        if isinstance(arg, omegaExpr):
            result = self.subs(omega / (2 * pi))
            return result.subs(arg, **assumptions)
        return super(fExpr, self).transform(arg, **assumptions)
        

class Yf(fExpr):

    """f-domain admittance"""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val, **assumptions):

        super(Yf, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Yt


class Zf(fExpr):

    """f-domain impedance"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, **assumptions):

        super(Zf, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Zt


class Hf(fExpr):

    """f-domain transfer function response."""

    quantity = 'Transfer function'
    units = ''

    def __init__(self, val, **assumptions):

        super(Hf, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Ht


class Vf(fExpr):

    """f-domain voltage (units V/Hz)."""

    quantity = 'Voltage spectrum'
    units = 'V/Hz'
    superkind = 'Voltage'        

    def __init__(self, val, **assumptions):

        super(Vf, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Vt


class If(fExpr):

    """f-domain current (units A/Hz)."""

    quantity = 'Current spectrum'
    units = 'A/Hz'
    superkind = 'Current'    

    def __init__(self, val, **assumptions):

        super(If, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = It


def fexpr(arg, **assumptions):
    """Create fExpr object.  If `arg` is fsym return f"""

    if arg is fsym:
        return f
    return fExpr(arg, **assumptions)
        
from .texpr import Ht, It, Vt, Yt, Zt, tExpr
f = fExpr('f')
