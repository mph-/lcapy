"""This module provides base transformer classes.

Copyright 2020--2022 Michael Hayes, UCECE

"""

from sympy import Eq, Symbol, Piecewise, S, Heaviside
from .state import state
from .sym import miscsymbol
from .utils import factor_const, remove_images
from .extrafunctions import UnitStep
from warnings import warn


class Transformer(object):

    name = 'undefined'
    is_inverse = False
    is_bilateral = False

    def __init__(self):

        self.cache = {}
        self.expr = None

    def clear_cache(self):

        self.cache = {}

    def error(self, message=''):
        if message == '':
            raise ValueError('Could not compute %s for %s' %
                             (self.name, self.expr))
        raise ValueError('Could not compute %s for %s: %s' %
                         (self.name, self.expr, message))

    def debug(self, message=''):
        if self._debug:
            print(self.name.capitalize() + ': ' + message)

    def simplify_term(self, expr, var):
        return expr

    def rewrite(self, expr, var):
        return expr

    def transform(self, expr, var, conjvar, evaluate=True, **kwargs):

        # Squirrel away original expression for error messages
        self.expr = expr
        self._debug = kwargs.pop('debug', False)

        if expr.is_Equality:
            return Eq(self.transform(expr.args[0], var, conjvar, evaluate, **kwargs),
                      self.transform(expr.args[1], var, conjvar, evaluate, **kwargs))

        # Handle Lcapy symbols for the variable and conjugate variable
        try:
            var = var.expr
        except:
            pass

        try:
            conjvar = conjvar.expr
        except:
            pass

        # The variable may have been created with different attributes,
        # say when using sym.sympify('DiracDelta(t)') since this will
        # default to assuming that t is complex.  So if the symbol has the
        # same representation, convert to the desired one.
        svar = Symbol(str(var))
        expr = expr.replace(svar, var)

        self.check(expr, var, conjvar, **kwargs)

        return self.doit(expr, var, conjvar, evaluate, **kwargs)

    def dummy_var(self, expr, dummy='m', level=0, **kwargs):
        """Create a dummy variable."""

        # Avoid i, j to reduce confusion with imaginary number
        # Avoid k, n since they are domain variables
        dummies = [dummy] + ['m', 'p', 'q', 'o']

        for dummy in dummies:
            if level == 0:
                nu = miscsymbol(dummy, **kwargs)
            else:
                nu = miscsymbol(dummy + '_%d' % level, **kwargs)
            if not expr.has(nu):
                return nu
        raise self.error(
            'Dummy variable conflict with symbols: %s' % ', '.join(dummies))


class BilateralForwardTransformer(Transformer):

    is_bilateral = True
    is_inverse = False

    def doit(self, expr, var, conjvar, evaluate=True, cache=True, **kwargs):

        const, expr = factor_const(expr, var)

        key = self.key(expr, var, conjvar, **kwargs)
        if cache and key in self.cache:
            return const * self.cache[key]

        expr = self.rewrite(expr, var)

        terms = expr.as_ordered_terms()
        result = 0

        for term in terms:
            sterm = self.simplify_term(term, var)
            ret = self.term(sterm, var, conjvar, **kwargs)
            result += ret

        self.cache[key] = result
        return const * result


class BilateralInverseTransformer(BilateralForwardTransformer):

    is_inverse = True


class UnilateralForwardTransformer(Transformer):

    is_bilateral = False

    def remove_heaviside(self, expr, var):

        rest = S.One
        for factor in expr.as_ordered_factors():
            if (factor.is_Function and
                factor.func in (Heaviside, UnitStep) and
                    factor.args[0] == var):
                pass
            else:
                rest *= factor
        return rest

    def simplify_term(self, expr, var):

        return self.remove_heaviside(expr, var)

    def doit(self, expr, var, conjvar, evaluate=True, **kwargs):

        # Unilateral transforms ignore expr for t < 0 so remove Piecewise.
        if expr.is_Piecewise and expr.args[0].args[1].has(var >= 0):
            expr = expr.args[0].args[0]

        if not evaluate:
            return self.noevaluate(expr, var, conjvar)

        const, expr = factor_const(expr, var)

        key = self.key(expr, var, conjvar, **kwargs)
        if key in self.cache:
            return const * self.cache[key]

        expr = self.rewrite(expr, var)

        terms = expr.as_ordered_terms()
        result = 0

        for term in terms:
            sterm = self.simplify_term(term, var)
            ret = self.term(sterm, var, conjvar, **kwargs)
            result += ret

        self.cache[key] = result
        return const * result


class UnilateralInverseTransformer(Transformer):

    is_bilateral = False
    is_inverse = True

    def make(self, var, const, cresult, uresult, **kwargs):

        result = const * (cresult + uresult)

        # The original signal can only be recovered from an inverse
        # unilateral Laplace transform if we have the region of
        # convergence.  In general, the result is unknown for t < 0
        # unless we know the result is causal.

        if not kwargs.get('causal', False):

            if uresult != 0:
                # Cannot determine result for var < 0
                result = Piecewise((result, var >= 0))

        return result

    def doit(self, expr, var, conjvar, evaluate=True, **kwargs):

        if not evaluate:
            return self.noevaluate(expr, var, conjvar)

        const, expr = factor_const(expr, var)

        key = self.key(expr, var, conjvar, **kwargs)
        if key in self.cache:
            return self.make(conjvar, const, *self.cache[key], **kwargs)

        expr = self.rewrite(expr, var)

        terms = expr.as_ordered_terms()

        uresult = 0
        cresult = 0

        for term in terms:
            sterm = self.simplify_term(term, var)
            cterm, uterm = self.term(sterm, var, conjvar, **kwargs)
            cresult += cterm
            uresult += uterm

        self.cache[key] = cresult, uresult
        return self.make(conjvar, const, *self.cache[key], **kwargs)
