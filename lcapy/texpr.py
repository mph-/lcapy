"""This module provides the tExpr class to represent time domain expressions.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from __future__ import division
from .expr import Expr
from .functions import exp
from .sym import fsym, ssym, tsym, j, oo, tausym
from .acdc import ACChecker, is_dc, is_ac, is_causal
from .laplace import laplace_transform
from .fourier import fourier_transform
from sympy import Heaviside, limit, Integral

__all__ = ('Ht', 'It', 'Vt', 'Yt', 'Zt')


class tExpr(Expr):

    """t-domain expression or symbol."""

    var = tsym
    domain_name = 'Time'
    domain_units = 's'

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)        
        assumptions['real'] = True
        super(tExpr, self).__init__(val, **assumptions)

        self._fourier_conjugate_class = fExpr
        self._laplace_conjugate_class = sExpr

        expr = self.expr        
        if check and expr.find(ssym) != set() and not expr.has(Integral):
            raise ValueError(
                't-domain expression %s cannot depend on s' % expr)
        if check and expr.find(fsym) != set() and not expr.has(Integral):
            raise ValueError(
                't-domain expression %s cannot depend on f' % expr)                            

    def infer_assumptions(self):

        self.assumptions['dc'] = False
        self.assumptions['ac'] = False
        self.assumptions['causal'] = False

        var = self.var
        if is_dc(self, var):
            self.assumptions['dc'] = True
            return

        if is_ac(self, var):
            self.assumptions['ac'] = True
            return

        if is_causal(self, var):
            self.assumptions['causal'] = True

    def laplace(self, evaluate=True, **assumptions):
        """Determine one-sided Laplace transform with 0- as the lower limit."""

        # The assumptions are required to help with the inverse Laplace
        # transform if required.
        self.infer_assumptions()

        assumptions = self.merge_assumptions(**assumptions)
        
        result = laplace_transform(self.expr, self.var, ssym, evaluate=evaluate)

        if hasattr(self, '_laplace_conjugate_class'):
            result = self._laplace_conjugate_class(result, **assumptions)
        else:
            result = sExpr(result, **assumptions)
        return result

    def LT(self, **assumptions):
        """Convert to s-domain.   This is an alias for laplace."""

        return self.laplace(**assumptions)
    
    def fourier(self, evaluate=True, **assumptions):
        """Attempt Fourier transform."""

        assumptions = self.merge_assumptions(**assumptions)
        
        result = fourier_transform(self.expr, self.var, fsym, evaluate=evaluate)

        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result, **assumptions)
        else:
            result = fExpr(result **self.assumptions)
        return result

    def FT(self, **assumptions):
        """Convert to f-domain.   This is an alias for fourier."""

        return self.fourier(**assumptions)    

    def phasor(self, **assumptions):

        check = ACChecker(self, t)
        if not check.is_ac:
            raise ValueError('Do not know how to convert %s to phasor' % self)
        phasor = Phasor(check.amp * exp(j * check.phase), omega=check.omega)
        return phasor

    def plot(self, t=None, **kwargs):
        """Plot the time waveform.  If t is not specified, it defaults to the
        range (-0.2, 2).  t can be a vector of specified instants, a
        tuple specifing the range, or a constant specifying the
        maximum value with the minimum value set to 0.

        kwargs include:
        axes - the plot axes to use otherwise a new figure is created
        xlabel - the x-axis label
        ylabel - the y-axis label
        xscale - the x-axis scaling, say for plotting as ms
        yscale - the y-axis scaling, say for plotting mV
        in addition to those supported by the matplotlib plot command.
        
        The plot axes are returned."""

        from .plot import plot_time
        return plot_time(self, t, **kwargs)

    def sample(self, t):

        """Return a discrete-time signal evaluated at time values specified by
        vector t. """

        return self.evaluate(t)

    def initial_value(self):
        """Determine value at t = 0. 
        See also pre_initial_value and post_initial_value"""

        return self.subs(0)

    def pre_initial_value(self):
        """Determine value at t = 0-.
        See also initial_value and post_initial_value"""

        return self.limit(self.var, 0, dir='-')

    def post_initial_value(self):
        """Determine value at t = 0+.
        See also pre_initial_value and initial_value"""

        return self.limit(self.var, 0, dir='+')    

    def final_value(self):
        """Determine value at t = oo."""

        return self.limit(self.var, oo)

    def remove_condition(self):
        """Remove the piecewise condition from the expression.
        See also force_causal."""

        if not self.is_conditional:
            return self
        expr = self.expr
        expr = expr.args[0].args[0]
        return self.__class__(expr)

    def force_causal(self):
        """Remove the piecewise condition from the expression
        and multiply by Heaviside function.  See also remove_condition."""

        if self.is_causal:
            return self
        
        expr = self.expr
        if self.is_conditional:
            expr = expr.args[0].args[0]            
        expr = expr * Heaviside(t)
        return self.__class__(expr)        

    def convolve(self, impulseresponse, commutate=False, **assumptions):
        """Convolve self with impulse response."""

        if not isinstance(impulseresponse, tExpr):
            raise ValueError('Expecting tExpr for impulse response')

        f1 = self.expr
        f2 = impulseresponse.expr
        if commutate:
            f1, f2 = f2, f1
        result = Integral(f1.subs(self.var, self.var - tausym) *
                          f2.subs(self.var, tausym),
                          (tausym, -oo, oo))
        return self.__class__(result, **assumptions)

    
class Yt(tExpr):

    """t-domain 'admittance' value."""

    units = 'siemens/s'

    def __init__(self, val, **assumptions):

        super(Yt, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Ys
        self._fourier_conjugate_class = Yf


class Zt(tExpr):

    """t-domain 'impedance' value."""

    units = 'ohms/s'

    def __init__(self, val, **assumptions):

        super(Zt, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Zs
        self._fourier_conjugate_class = Zf


class Vt(tExpr):

    """t-domain voltage (units V)."""

    quantity = 'Voltage'
    units = 'V'
    superkind = 'Voltage'        

    def __init__(self, val, **assumptions):

        super(Vt, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Vs
        self._fourier_conjugate_class = Vf

    def cpt(self):
        from .oneport import V
        return V(self)

    
class It(tExpr):

    """t-domain current (units A)."""

    quantity = 'Current'
    units = 'A'
    superkind = 'Current'    

    def __init__(self, val, **assumptions):

        super(It, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Is
        self._fourier_conjugate_class = If

    def cpt(self):
        from .oneport import I
        return I(self)
        

class Ht(tExpr):

    """impulse response"""

    quantity = 'Impulse response'
    units = '1/s'

    def __init__(self, val, **assumptions):

        super(Ht, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Hs
        self._fourier_conjugate_class = Hf

def texpr(arg, **assumptions):
    """Create tExpr object.  If `arg` is tsym return t"""

    if arg is tsym:
        return t
    return tExpr(arg, **assumptions)

from .sexpr import Hs, Is, Vs, Ys, Zs, sExpr
from .fexpr import Hf, If, Vf, Yf, Zf, fExpr
from .phasor import Phasor
t = tExpr('t')
