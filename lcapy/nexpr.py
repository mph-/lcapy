"""This module provides the DiscreteTimeDomainExpression class to
represent discrete-time expressions.

Copyright 2020--2021 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import DiscreteTimeDomain
from .seqexpr import SequenceExpression
from .sequence import Sequence
from .functions import exp
from .sym import j, oo, pi, fsym
from .dsym import nsym, ksym, zsym, dt
from .ztransform import ztransform
from .dft import DFT
from sympy import Sum, summation, limit


__all__ = ('nexpr', )

class DiscreteTimeDomainExpression(DiscreteTimeDomain, SequenceExpression):
    """Discrete time expression or symbol."""

    var = nsym


    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)
        assumptions['real'] = True        
        
        super(DiscreteTimeDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr
        if check and expr.find(zsym) != set() and not expr.has(Sum):
            raise ValueError(
                'n-domain expression %s cannot depend on z' % expr)
        if check and expr.find(ksym) != set() and not expr.has(Sum):
            raise ValueError(
                'n-domain expression %s cannot depend on k' % expr)            

    def _mul_compatible_domains(self, x):

        if self.domain == x.domain:
            return True        

        return x.is_constant_domain

    def _div_compatible_domains(self, x):

        if self.domain == x.domain:
            return True
        
        return x.is_constant_domain

    def as_expr(self):
        return DiscreteTimeDomainExpression(self)

    def differentiate(self):
        """First order difference."""

        result = (self.expr - self.subs(n - 1).expr) / dt
        return self.__class__(result, **self.assumptions)

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

        assumptions = self.assumptions.merge_and_infer(self, **assumptions)
        result = ztransform(self.expr, self.var, zsym, evaluate)
        return self.change(result, domain='Z', **assumptions)

    def ZT(self, **assumptions):
        return self.ztransform(**assumptions)

    def plot(self, ni=None, **kwargs):
        """Plot the sequence.  If `ni` is not specified, it defaults to the
        range (-20, 20).  `ni` can be a vector of specified sequence
        indices, a tuple specifing the range, or a constant specifying
        the maximum value with the minimum value set to 0.

        kwargs include:
        axes - the plot axes to use otherwise a new figure is created
        xlabel - the x-axis label
        ylabel - the y-axis label
        xscale - the x-axis scaling, say for plotting as ms
        yscale - the y-axis scaling, say for plotting mV
        in addition to those supported by the matplotlib plot command.
        
        The plot axes are returned.

        """

        if ni is None:
            ni = (-20, 20)
        
        from .plot import plot_sequence
        return plot_sequence(self, ni, **kwargs)

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
        return self.change(result, domain='discrete fourier')
    
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

        assumptions = self.assumptions.merge_and_infer(self, **assumptions)
        if assumptions.is_causal:
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

    def remove_condition(self):
        """Remove the piecewise condition from the expression."""

        if not self.is_conditional:
            return self
        expr = self.expr
        expr = expr.args[0].args[0]
        return self.__class__(expr)
    
    
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

    return DiscreteTimeDomainExpression(arg, **assumptions)


from .expressionclasses import expressionclasses

classs = expressionclasses.register('discrete time', DiscreteTimeDomainExpression)

n = DiscreteTimeDomainExpression('n')
