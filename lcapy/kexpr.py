"""This file provides the kExpr class to represent k-domain (discrete Fourier
domain) expressions.

Copyright 2020 Michael Hayes, UCECE

"""

from __future__ import division
from .fourier import inverse_fourier_transform
from .expr import Expr
from .dsym import nsym, ksym, zsym, dt, df


class kExpr(Expr):

    """Fourier domain expression or symbol."""

    var = ksym
    domain_name = 'Frequency'
    domain_units = 'Hz'

    def __init__(self, val, **assumptions):

        assumptions['real'] = True
        super(kExpr, self).__init__(val, **assumptions)
        # Define when class defined.
        self._fourier_conjugate_class = nexpr

        if self.expr.find(zsym) != set():
            raise ValueError(
                'k-domain expression %s cannot depend on z' % self.expr)
        if self.expr.find(nsym) != set():
            raise ValueError(
                'f-domain expression %s cannot depend on t' % self.expr)

    def inverse_fourier(self):
        """Attempt inverse Fourier transform."""

        result = inverse_fourier_transform(self.expr, self.var, nsym)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        else:
            result = nexpr(result)
        return result

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


class Yf(kExpr):

    """f-domain admittance"""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val, **assumptions):

        super(Yf, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Yn


class Zf(kExpr):

    """f-domain impedance"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, **assumptions):

        super(Zf, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Zn


class Hf(kExpr):

    """f-domain transfer function response."""

    quantity = 'Transfer function'
    units = ''

    def __init__(self, val, **assumptions):

        super(Hf, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Hn


class Vf(kExpr):

    """f-domain voltage (units V/Hz)."""

    quantity = 'Voltage spectrum'
    units = 'V/Hz'

    def __init__(self, val, **assumptions):

        super(Vf, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = Vn


class If(kExpr):

    """f-domain current (units A/Hz)."""

    quantity = 'Current spectrum'
    units = 'A/Hz'

    def __init__(self, val, **assumptions):

        super(If, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = In


def kexpr(arg):
    """Create kExpr object.  If `arg` is ksym return k"""

    if arg is ksym:
        return k
    return kExpr(arg)
        
from .nexpr import Hn, In, Vn, Yn, Zn, nexpr
k = kExpr('k')
