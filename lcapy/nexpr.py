"""This module provides the nExpr class to represent discrete-time expressions.

Copyright 2020 Michael Hayes, UCECE

"""

from __future__ import division
from .seqexpr import seqExpr
from .functions import exp
from .sym import j, oo, pi, fsym
from .dsym import nsym, ksym, zsym, dt
from .acdc import ACChecker, is_dc, is_ac, is_causal
from .ztransform import ztransform
from .dft import DFT
from sympy import Sum, summation, limit


__all__ = ('Hn', 'In', 'Vn', 'Yn', 'Zn')


class nExpr(seqExpr):

    """t-domain expression or symbol."""

    var = nsym
    domain_name = 'Sample'
    domain_units = ''

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)
        
        super(nExpr, self).__init__(val, **assumptions)

        self._discrete_fourier_conjugate_class = kExpr
        self._ztransform_conjugate_class = zExpr

        expr = self.expr
        if check and expr.find(zsym) != set() and not expr.has(Sum):
            raise ValueError(
                'n-domain expression %s cannot depend on z' % expr)
        if check and expr.find(ksym) != set() and not expr.has(Sum):
            raise ValueError(
                'n-domain expression %s cannot depend on k' % expr)            

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

        # TODO, get SymPy to optimize this case.
        expr = self.expr
        const, expr = factor_const(expr, nsym)
        if expr.is_Function and expr.func == UnitImpulse:
            return dt * u(expr.args[0]) * const
        
        msym = sympify('m', real=True)
        result = dt * summation(self.subs(msym).expr, (msym, -oo, nsym))
        
        return self.__class__(result, **self.assumptions)

    def ztransform(self, evaluate=True, **assumptions):
        """Determine one-sided z-transform."""

        self.infer_assumptions()

        assumptions = self.merge_assumptions(**assumptions)

        result = ztransform(self.expr, self.var, zsym, evaluate)

        if hasattr(self, '_ztransform_conjugate_class'):
            result = self._ztransform_conjugate_class(result, **assumptions)
        else:
            result = zExpr(result, **assumptions)
        return result

    def ZT(self, **assumptions):
        return self.ztransform(**assumptions)

    def plot(self, n=None, **kwargs):
        """Plot the sequence.  If n is not specified, it defaults to the
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

        return self.__class__(limit(self.expr, self.var, oo))

    def DFT(self, N=None, evaluate=True):

        from .kexpr import k

        if N is None:
            from .sym import sympify
            
            N = sympify('N')

        result = DFT(self.expr, nsym, ksym, N, evaluate=evaluate)

        if hasattr(self, '_discrete_fourier_conjugate_class'):
            result = self._discrete_fourier_conjugate_class(result)
        else:
            result = kExpr(result, check=False)
            
        return result
    
    def delay(self,m):
        """Delay signal by m samples."""

        return self.subs(n - m)

    def extent(self, n1=-100, n2=100):
        """Determine extent of the signal.

        For example, nexpr([1, 1]).extent() = 2
                     nexpr([1, 0, 1]).extent() = 3
                     nexpr([0, 1, 0, 1]).extent() = 3

        This performs a search between n=n1 and n=n2."""

        return self.seq((n1, n2)).extent()

    def discrete_time_fourier_transform(self, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""

        self.infer_assumptions()

        assumptions = self.merge_assumptions(**assumptions)
        
        if assumptions.get('causal', False):        
            return self.ZT(**assumptions).DTFT()

        from .fexpr import fexpr

        # TODO, handle evaluation keyword
        return fexpr(summation(self.subs(nsym).expr * exp(-2 * j * pi * fsym * nsym * dt), (nsym, -oo, oo)))

    def DTFT(self, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""
    
        return self.discrete_time_fourier_transform(**assumptions)

    def difference_equation(self, input='x', output='y', form='iir'):
        """Create difference equation from impulse response.

        form can be 'fir' or 'iir' ('direct form I').
        """

        H = self.ZT()
        return H.difference_equation(input, output, form)

    
class Yn(nExpr):

    """t-domain 'admittance' value."""

    units = 'siemens/s'

    def __init__(self, val, **assumptions):

        super(Yn, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Yz
        self._discrete_fourier_conjugate_class = Yk


class Zn(nExpr):

    """t-domain 'impedance' value."""

    units = 'ohms/s'

    def __init__(self, val, **assumptions):

        super(Zn, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Zz
        self._discrete_fourier_conjugate_class = Zk


class Vn(nExpr):

    """t-domain voltage (units V)."""

    quantity = 'Voltage'
    units = 'V'

    def __init__(self, val, **assumptions):

        super(Vn, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Vz
        self._discrete_fourier_conjugate_class = Vk

class In(nExpr):

    """t-domain current (units A)."""

    quantity = 'Current'
    units = 'A'

    def __init__(self, val, **assumptions):

        super(In, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Iz
        self._discrete_fourier_conjugate_class = Ik


class Hn(nExpr):

    """impulse response"""

    quantity = 'Impulse response'
    units = ''

    def __init__(self, val, **assumptions):

        super(Hn, self).__init__(val, **assumptions)
        self._ztransform_conjugate_class = Hz
        self._discrete_fourier_conjugate_class = Hk

def nexpr(arg, **assumptions):
    """Create nExpr object.  If `arg` is nsym return n"""

    if arg is nsym:
        return n

    from .seq import seq
    
    if isinstance(arg, str) and arg.startswith('{'):
        return seq(arg)
    
    from numpy import ndarray
    
    if isinstance(arg, (list, ndarray)):
        return Sequence(arg, var=n).as_impulses()

    return nExpr(arg, **assumptions)

from .zexpr import Hz, Iz, Vz, Yz, Zz, zExpr
from .kexpr import Hk, Ik, Vk, Yk, Zk, kExpr
n = nExpr('n')
