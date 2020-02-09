"""This file provides the nExpr class to represent discrete-time expressions.

Copyright 2020 Michael Hayes, UCECE

"""

from __future__ import division
from .expr import Expr
from .functions import exp
from .sym import j, oo, pi
from .dsym import nsym, ksym, zsym, dt
from .acdc import ACChecker, is_dc, is_ac, is_causal
from .ztransform import ztransform
from .sequence import Sequence
from sympy import Sum, summation


__all__ = ('Hn', 'In', 'Vn', 'Yn', 'Zn')


class nExpr(Expr):

    """t-domain expression or symbol."""

    var = nsym
    domain_name = 'Sample'
    domain_units = ''

    def __init__(self, val, **assumptions):

        assumptions['real'] = True
        super(nExpr, self).__init__(val, **assumptions)

        #self._fourier_conjugate_class = fExpr
        self._ztransform_conjugate_class = zExpr

        if self.expr.find(zsym) != set():
            raise ValueError(
                'n-domain expression %s cannot depend on s' % self.expr)

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

    def merge_assumptions(self, **assumptions):
        
        new_assumptions = self.assumptions.copy()
        new_assumptions.update(assumptions)
        return new_assumptions

    def differentiate(self):
        """First order difference."""

        return self.__class__((self.expr - self.subs(n - 1).expr) / dt, **self.assumptions)

    def integrate(self):
        """First order integration."""

        from .sym import sympify
        from .utils import factor_const
        from .functions import UnitImpulse, u

        expr = self.expr
        const, expr = factor_const(expr, nsym)
        if expr.is_Function and expr.func == UnitImpulse:
            return dt * u(expr.args[0]) * const
        
        msym = sympify('m', real=True)
        result = dt * summation(self.subs(msym).expr, (msym, -oo, nsym))
        
        return self.__class__(result, **self.assumptions)

    def ztransform(self, **assumptions):
        """Determine one-sided z-transform."""

        self.infer_assumptions()

        assumptions = self.merge_assumptions(**assumptions)

        result = ztransform(self.expr, self.var, zsym)

        if hasattr(self, '_ztransform_conjugate_class'):
            result = self._ztransform_conjugate_class(result, **assumptions)
        else:
            result = zExpr(result, **assumptions)
        return result

    def ZT(self, **assumptions):
        return self.ztransform(**assumptions)

    def plot(self, n=None, **kwargs):
        """Plot the time waveform.  If n is not specified, it defaults to the
        range (-20, 20).  n can be a vector of specified instants, a
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

        if n is None:
            n = (-20, 20)
        
        from .plot import plot_sequence
        return plot_sequence(self, n, **kwargs)

    def initial_value(self):
        """Determine value at n = 0."""

        return self.subs(0)

    def final_value(self):
        """Determine value at n = oo."""

        return self.__class__(sym.limit(self.expr, self.var, oo))

    def seq(self, nvals=None, evaluate=False):

        if nvals is None:
            nvals = (-10, 10)
        if isinstance(nvals, tuple):
            nvals = range(nvals[0], nvals[1] + 1)

        nvals = list(nvals)
        v = self(nvals)

        return Sequence(v, nvals, evaluate, self.var)

    def DFT(self, N=None, evaluate=True):

        from .kexpr import k

        if N is None:
            from .expr import expr
            from .sym import sympify
            
            N = sympify('N')

        foo = self.expr * exp(-2 * j * pi * nsym * ksym / N)

        if evaluate:
            result = summation(foo, (nsym, 0, N))                        
        else:
            result = Sum(foo, (nsym, 0, N))

        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        else:
            result = kExpr(result, check=False)
            
        return result
    
    
class Yn(nExpr):

    """t-domain 'admittance' value."""

    units = 'siemens/s'

    def __init__(self, val, **assumptions):

        super(Yn, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Yz
        self._fourier_conjugate_class = Yk


class Zn(nExpr):

    """t-domain 'impedance' value."""

    units = 'ohms/s'

    def __init__(self, val, **assumptions):

        super(Zn, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Zz
        self._fourier_conjugate_class = Zk


class Vn(nExpr):

    """t-domain voltage (units V)."""

    quantity = 'Voltage'
    units = 'V'

    def __init__(self, val, **assumptions):

        super(Vn, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Vz
        self._fourier_conjugate_class = Vk

class In(nExpr):

    """t-domain current (units A)."""

    quantity = 'Current'
    units = 'A'

    def __init__(self, val, **assumptions):

        super(In, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Iz
        self._fourier_conjugate_class = Ik


class Hn(nExpr):

    """impulse response"""

    quantity = 'Impulse response'
    units = ''

    def __init__(self, val, **assumptions):

        super(Hn, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Hz
        self._fourier_conjugate_class = Hk

def nexpr(arg):
    """Create nExpr object.  If `arg` is nsym return n"""

    if arg is nsym:
        return n
    return nExpr(arg)

from .zexpr import Hz, Iz, Vz, Yz, Zz, zExpr
from .kexpr import Hk, Ik, Vk, Yk, Zk, kExpr
n = nExpr('n')
