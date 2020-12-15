"""This module provides the FourierDomainExpression class to represent f-domain (Fourier
domain) expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .fourier import inverse_fourier_transform
from .expr import Expr, expr
from .sym import fsym, ssym, tsym, pi
from sympy import Integral, Expr as symExpr

class FourierDomainExpression(Expr):

    """Fourier domain expression or symbol."""

    var = fsym
    domain_name = 'Frequency'
    domain_units = 'Hz'

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)        
        assumptions['real'] = True
        super(FourierDomainExpression, self).__init__(val, **assumptions)
        # Define when class defined.
        self._fourier_conjugate_class = TimeDomainExpression

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
            result = TimeDomainExpression(result)
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

        from .omegaexpr import AngularFourierDomainExpression, omega

        arg = expr(arg)        
        if isinstance(arg, AngularFourierDomainExpression):
            result = self.subs(omega / (2 * pi))
            return result.subs(arg, **assumptions)
        return super(FourierDomainExpression, self).transform(arg, **assumptions)
        

class FourierDomainAdmittance(FourierDomainExpression):

    """f-domain admittance"""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val, **assumptions):

        super(FourierDomainAdmittance, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainAdmittance


class FourierDomainImpedance(FourierDomainExpression):

    """f-domain impedance"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, **assumptions):

        super(FourierDomainImpedance, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainImpedance


class FourierDomainTransferFunction(FourierDomainExpression):

    """f-domain transfer function response."""

    quantity = 'Transfer function'
    units = ''

    def __init__(self, val, **assumptions):

        super(FourierDomainTransferFunction, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainImpulseResponse


class FourierDomainVoltage(FourierDomainExpression):

    """f-domain voltage (units V/Hz)."""

    quantity = 'Voltage spectrum'
    units = 'V/Hz'

    def __init__(self, val, **assumptions):

        super(FourierDomainVoltage, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainVoltage

    def __mul__(self, x):
        """Multiply"""

        if isinstance(x, FourierDomainAdmittance):
            return FourierDomainCurrent(super(FourierDomainVoltage, self).__mul__(x))
        if isinstance(x, (ConstantExpression, FourierDomainExpression, symExpr, int, float, complex)):
            return super(FourierDomainVoltage, self).__mul__(x)
        self._incompatible(x, '*')

    def __truediv__(self, x):
        """Divide"""

        if isinstance(x, FourierDomainImpedance):
            return FourierDomainCurrent(super(FourierDomainVoltage, self).__truediv__(x))
        if isinstance(x, FourierDomainVoltage):
            return FourierDomainTransferFunction(super(FourierDomainVoltage, self).__truediv__(x))                
        if isinstance(x, (ConstantExpression, FourierDomainExpression, symExpr, int, float, complex)):
            return super(FourierDomainVoltage, self).__truediv__(x)
        self._incompatible(x, '/')                

        
class FourierDomainCurrent(FourierDomainExpression):

    """f-domain current (units A/Hz)."""

    quantity = 'Current spectrum'
    units = 'A/Hz'

    def __init__(self, val, **assumptions):

        super(FourierDomainCurrent, self).__init__(val, **assumptions)
        self._fourier_conjugate_class = TimeDomainCurrent

    def __mul__(self, x):
        """Multiply"""

        if isinstance(x, FourierDomainImpedance):
            return FourierDomainVoltage(super(FourierDomainCurrent, self).__mul__(x))            
        if isinstance(x, (ConstantExpression, FourierDomainExpression, symExpr, int, float, complex)):
            return super(FourierDomainCurrent, self).__mul__(x)
        self._incompatible(x, '*')        

    def __truediv__(self, x):
        """Divide"""

        if isinstance(x, FourierDomainAdmittance):
            return FourierDomainVoltage(super(FourierDomainCurrent, self).__truediv__(x))
        if isinstance(x, FourierDomainCurrent):
            return FourierDomainTransferFunction(super(FourierDomainCurrent, self).__truediv__(x))                
        if isinstance(x, (ConstantExpression, FourierDomainExpression, symExpr, int, float, complex)):
            return super(FourierDomainCurrent, self).__truediv__(x)
        self._incompatible(x, '/')                        

        
def fexpr(arg, **assumptions):
    """Create FourierDomainExpression object.  If `arg` is fsym return f"""

    if arg is fsym:
        return f
    return FourierDomainExpression(arg, **assumptions)
        
from .texpr import TimeDomainImpulseResponse, TimeDomainCurrent, TimeDomainVoltage, TimeDomainAdmittance, TimeDomainImpedance, TimeDomainExpression
from .cexpr import ConstantExpression
f = FourierDomainExpression('f')
