"""This module provides the DiscreteTimeDomainExpression class to
represent discrete-time expressions.

Copyright 2020--2022 Michael Hayes, UCECE

"""

from __future__ import division
from .domains import DiscreteTimeDomain
from .functions import UnitStep
from .sym import oo, fsym
from .sym import nsym, ksym, zsym, dt
from .ztransform import ztransform
from .dft import DFT
from .seqexpr import SequenceExpression
from .nseq import DiscreteTimeDomainSequence, nseq
from sympy import Sum, summation, limit


__all__ = ('nexpr', )


class DiscreteTimeDomainExpression(DiscreteTimeDomain, SequenceExpression):
    """Discrete-time expression or symbol."""

    var = nsym
    seqcls = DiscreteTimeDomainSequence

    def __init__(self, val, **assumptions):

        check = assumptions.pop('check', True)
        if 'integer' not in assumptions:
            assumptions['real'] = True

        super(DiscreteTimeDomainExpression, self).__init__(val, **assumptions)

        expr = self.expr
        if check and expr.has(zsym) and not expr.has(Sum):
            raise ValueError(
                'n-domain expression %s cannot depend on z' % expr)
        if check and expr.has(ksym) and not expr.has(Sum):
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

        from .sym import miscsymbol
        from .utils import factor_const
        from .extrafunctions import UnitImpulse
        from .functions import u

        # TODO, get SymPy to optimize this case.
        expr = self.expr
        const, expr = factor_const(expr, nsym)
        if expr.is_Function and expr.func == UnitImpulse:
            return dt * u(expr.args[0]) * const

        msym = miscsymbol('m', integer=True)
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

    def DFT(self, N=None, evaluate=True, piecewise=False):
        """Determine DFT.

        `N` needs to be a positive integer symbol or a str specifying
        the extent of the DFT.  By default `N` is defined as 'N'."""

        from .sym import miscsymbol

        if N is None:
            N = miscsymbol('N', integer=True, positive=True)
        elif isinstance(N, str):
            N = miscsymbol(N, integer=True, positive=True)

        result = DFT(self.expr, nsym, ksym, N, evaluate=evaluate,
                     piecewise=piecewise)
        result = self.change(result, domain='discrete fourier')
        result = result.simplify_unit_impulse()
        return result

    def delay(self, m):
        """Delay signal by m samples."""

        return self.subs(n - m)

    def extent(self, n1=-100, n2=100):
        """Determine extent of the signal.

        For example, nexpr([1, 1]).extent() = 2
                     nexpr([1, 0, 1]).extent() = 3
                     nexpr([0, 1, 0, 1]).extent() = 3

        This performs a search between n=n1 and n=n2."""

        return self.seq((n1, n2)).extent()

    def discrete_time_fourier_transform(self, var=None,
                                        images=oo, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform.

        Use `images = 0` to avoid the infinite number of spectral images.
        """

        return self.DTFT(var, images, **assumptions)

    def DTFT(self, var=None, images=oo, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform.

        By default this returns the DTFT in terms of `f`.  Use
        `.DTFT(w)` to get the angular frequency form, `.DTFT(F)` to
        get the normalised frequency form, or `.DTFT(W)` to get the
        normalised angular frequency form.

        Use `images = 0` to avoid the infinite number of spectral images.

        """

        from .symbols import f, omega, Omega, F
        from .fexpr import fexpr
        from .dtft import DTFT

        if var is None:
            var = f
        if id(var) not in (id(f), id(F), id(omega), id(Omega)):
            raise ValueError(
                'DTFT requires var to be f, F, omega, or Omega`, not %s' % var)

        dtft = DTFT(self.expr, self.var, fsym, images=images)

        result = fexpr(dtft)(var)
        result = result.simplify_dirac_delta()
        result = result.simplify_heaviside()
        result = result.simplify_rect()

        # There is a bug in SymPy when simplifying Sum('X(n - m)', (m, -oo, oo))
        # result = result.simplify()
        result = result.cancel_terms()

        return result

    def norm_angular_fourier(self, **assumptions):

        from .normomegaexpr import Omega
        return self.DTFT()(Omega)

    def difference_equation(self, inputsym='x', outputsym='y', form='iir'):
        """Create difference equation from impulse response.

        `form` can be 'fir' or 'iir' ('direct form I').
        """

        H = self.ZT()
        return H.difference_equation(inputsym, outputsym, form)

    def remove_condition(self):
        """Remove the piecewise condition from the expression."""

        if not self.is_conditional:
            return self
        expr = self.expr
        expr = expr.args[0].args[0]
        return self.__class__(expr)

    def force_causal(self):
        """Remove the piecewise condition from the expression
        and multiply by unit-step function.  See also remove_condition."""

        if self.is_causal:
            return self

        expr = self.expr
        if self.is_conditional:
            expr = expr.args[0].args[0]
        expr = expr * UnitStep(n.var)
        return self.__class__(expr)

    def zdomain(self, **assumptions):
        return self.ZT(**assumptions)

    def discrete_frequency(self, **assumptions):
        return self.DFT(**assumptions)

    def discrete_time(self, **assumptions):
        return self

    def fourier(self, **assumptions):
        return self.DTFT(**assumptions)

    def angular_fourier(self, **assumptions):
        from .symbols import omega

        return self.DTFT(omega, **assumptions)

    def norm_fourier(self, **assumptions):
        from .symbols import F

        return self.DTFT(F, **assumptions)


def nexpr(arg, **assumptions):
    """Create nExpr object.  If `arg` is nsym return n"""

    from .expr import Expr

    if arg is nsym:
        return n

    if isinstance(arg, Expr):
        if assumptions == {}:
            return arg
        return arg.__class__(arg, **assumptions)

    if isinstance(arg, str) and arg.startswith('{'):
        return nseq(arg)

    from numpy import ndarray

    if isinstance(arg, (list, ndarray)):
        return DiscreteTimeDomainSequence(arg, var=n).as_impulses()

    return DiscreteTimeDomainExpression(arg, **assumptions)


from .expressionclasses import expressionclasses  # nopep8

classes = expressionclasses.register(
    'discrete time', DiscreteTimeDomainExpression)
DiscreteTimeDomainVoltage = classes['voltage']
DiscreteTimeDomainCurrent = classes['current']

n = DiscreteTimeDomainExpression('n', integer=True)
