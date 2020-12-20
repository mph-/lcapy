"""This module provides the FourierDomainExpression class to represent f-domain (Fourier
domain) expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .fourier import inverse_fourier_transform
from .expr import Expr, expr
from .sym import fsym, ssym, tsym, pi
from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
from .admittancemixin import AdmittanceMixin
from .impedancemixin import ImpedanceMixin
from .transfermixin import TransferMixin
from sympy import Integral, Expr as symExpr

class FourierDomainExpression(Expr):

    """Fourier domain expression or symbol."""

    var = fsym
    domain_name = 'Frequency'
    domain_units = 'Hz'
    is_fourier_domain = True    

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

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        from .symbols import omega
        
        if assumptions.get('causal', self.is_causal):
            result = self.subs(omega / (2 * pi))
        else:
            result = self.time(**assumptions).angular_fourier(**assumptions)

        return self.wrap(result)            

    def inverse_fourier(self, evaluate=True, **assumptions):
        """Attempt inverse Fourier transform."""

        result = inverse_fourier_transform(self.expr, self.var, tsym, evaluate=evaluate)

        return self.wrap(TimeDomainExpression(result, **assumptions))        

    def IFT(self, **assumptions):
        """Convert to t-domain.   This is an alias for inverse_fourier."""

        return self.inverse_fourier(**assumptions)    
    
    def time(self, **assumptions):
        return self.inverse_fourier(**assumptions)

    def angular_fourier(self, **assumptions):
        """Convert to angular Fourier domain."""
        from .symbols import omega
        
        result = self.subs(omega / (2 * pi))
        return self.wrap(result)            
    
    def laplace(self, **assumptions):
        """Determine one-side Laplace transform with 0- as the lower limit."""

        from .symbols import s, j

        if self.is_causal:
            result = self.subs(s / (j * 2 * pi))
        else:
            result = self.time.laplace()

        return self.wrap(LaplaceDomainExpression(result, **assumptions))
    
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


class FourierDomainAdmittance(AdmittanceMixin, FourierDomainExpression):
    """f-domain admittance"""

    quantity = 'Admittance'
    units = 'siemens'


class FourierDomainImpedance(ImpedanceMixin, FourierDomainExpression):
    """f-domain impedance"""

    quantity = 'Impedance'
    units = 'ohms'


class FourierDomainTransferFunction(TransferMixin, FourierDomainExpression):
    """f-domain transfer function response."""

    quantity = 'Transfer function'
    units = ''


class FourierDomainVoltage(VoltageMixin, FourierDomainExpression):
    """f-domain voltage (units V/Hz)."""

    quantity = 'Voltage spectrum'
    units = 'V/Hz'

        
class FourierDomainCurrent(CurrentMixin, FourierDomainExpression):
    """f-domain current (units A/Hz)."""

    quantity = 'Current spectrum'
    units = 'A/Hz'

        
def fexpr(arg, **assumptions):
    """Create FourierDomainExpression object.  If `arg` is fsym return f"""

    if arg is fsym:
        return f
    return FourierDomainExpression(arg, **assumptions)
        
from .texpr import TimeDomainExpression
from .sexpr import LaplaceDomainExpression
from .omegaexpr import AngularFourierDomainExpression
from .cexpr import ConstantExpression
f = FourierDomainExpression('f')
