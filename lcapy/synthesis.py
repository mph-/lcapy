"""
This module implements experimental network synthesis.

Copyright 2020 Michael Hayes, UCECE
"""

from .oneport import L, C, R, G, parallel, series
from .sexpr import s, LaplaceDomainImpedance, LaplaceDomainExpression
from .expr import expr

# Should check that args to L, C, R, G contains s and raise
# exception since the circuit is not realisable.

# Could make this a mixin for Zs and Ys but would need a Z flavour
# and a Y flavour.


class Synthesis(object):

    def seriesRL(self, lexpr):
        """Z = s * L + R"""

        if lexpr == 0:
            # Perhaps return Wire object
            return None

        var = lexpr.var
        lexpr = lexpr.partfrac()
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = R(a)

        a = d.pop(var, None)
        if a is not None:
            net = series(net, L(a))

        if d != {}:
            raise ValueError('Not series RL')

        return net

    def seriesRC(self, lexpr):
        """Z = R + 1 / (s * C)"""

        if lexpr == 0:
            raise ValueError('Not series RC')

        var = lexpr.var
        lexpr = lexpr.partfrac()
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = R(a)

        a = d.pop(1 / var, None)
        if a is not None:
            net = series(net, C(1 / a))

        if d != {}:
            raise ValueError('Not series RC')

        return net

    def seriesGC(self, lexpr):
        """Z = 1 / G + 1 / (s * C)"""

        if lexpr == 0:
            raise ValueError('Not series GC')

        var = lexpr.var
        lexpr = lexpr.partfrac()
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = G(1 / a)

        a = d.pop(1 / var, None)
        if a is not None:
            net = series(net, C(1 / a))

        if d != {}:
            raise ValueError('Not series GC')

        return net

    def seriesLC(self, lexpr):
        """Z = s * L + 1 / (s * C)"""

        if lexpr == 0:
            raise ValueError('Not series LC')

        var = lexpr.var
        lexpr = lexpr.partfrac()
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, 0)
        if a != 0:
            raise ValueError('Not series LC')

        a = d.pop(1 / var, None)
        if a is not None:
            net = C(1 / a)

        a = d.pop(var, None)
        if a is not None:
            net = series(net, L(a))

        if d != {}:
            raise ValueError('Not series LC')

        return net

    def seriesRLC(self, lexpr):
        """Z = s * L + R + 1 / (s * C)"""

        if lexpr == 0:
            raise ValueError('Not series RLC')

        var = lexpr.var
        lexpr = lexpr.partfrac()
        d = lexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = R(a)

        a = d.pop(1 / var, None)
        if a is not None:
            net = series(net, C(1 / a))

        a = d.pop(var, None)
        if a is not None:
            net = series(net, L(a))

        if d != {}:
            raise ValueError('Not series RLC')

        return net

    def parallelRL(self, lexpr):
        """Y = 1 / (s * L) + 1 / R"""

        if lexpr == 0:
            # Perhaps return Wire object
            return None

        var = lexpr.var
        yexpr = (1 / lexpr).partfrac()
        d = yexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = R(1 / a)

        a = d.pop(1 / var, None)
        if a is not None:
            net = parallel(net, L(1 / a))

        if d != {}:
            raise ValueError('Not parallel RL')

        return net

    def parallelRC(self, lexpr):
        """Y = s * C + 1 / R"""

        if lexpr == 0:
            raise ValueError('Not parallel RC')

        var = lexpr.var
        yexpr = (1 / lexpr).partfrac()
        d = yexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = R(1 / a)

        a = d.pop(var, None)
        if a is not None:
            net = parallel(net, C(a))

        if d != {}:
            raise ValueError('Not parallel RC')

        return net

    def parallelGC(self, lexpr):
        """Y = s * C + G"""

        if lexpr == 0:
            raise ValueError('Not parallel GC')

        var = lexpr.var
        yexpr = (1 / lexpr).partfrac()
        d = yexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, None)
        if a is not None:
            net = G(a)

        a = d.pop(var, None)
        if a is not None:
            net = parallel(net, C(a))

        if d != {}:
            raise ValueError('Not parallel GC')

        return net

    def parallelLC(self, lexpr):
        """Y = s * C + 1 / (s * L)"""

        if lexpr == 0:
            raise ValueError('Not parallel LC')

        var = lexpr.var
        yexpr = (1 / lexpr).partfrac()
        d = yexpr.expr.collect(var, evaluate=False)

        net = None
        a = d.pop(1, 0)
        if a != 0:
            raise ValueError('Not parallel LC')

        a = d.pop(var, None)
        if a is not None:
            net = C(a)

        a = d.pop(1 / var, None)
        if a is not None:
            net = parallel(net, L(1 / a))

        if d != {}:
            raise ValueError('Not parallel LC')

        return net

    def parallelRLC(self, lexpr):
        """Y = s * C + 1 / R + 1 / (s * L)"""

        if lexpr == 0:
            raise ValueError('Not parallel RLC')

        var = lexpr.var
        yexpr = (1 / lexpr).partfrac()
        d = yexpr.expr.collect(var, evaluate=False)

        Rnet = None
        Lnet = None
        Cnet = None
        a = d.pop(1, None)
        if a is not None:
            Rnet = R(1 / a)

        a = d.pop(var, None)
        if a is not None:
            Cnet = C(a)

        a = d.pop(1 / var, None)
        if a is not None:
            Lnet = L(1 / a)

        if d != {}:
            raise ValueError('Not parallel RLC')

        return parallel(Rnet, Lnet, Cnet)

    def RLC(self, lexpr):
        """Handle series RLC or parallel RLC.  This will fail on R | L + C."""

        try:
            return self.seriesRLC(lexpr)
        except:
            return self.parallelRLC(lexpr)

    def parallelseriesRL(self):
        """(R('R1') + L('L1')) | (R('R2') + L('L2'))"""
        # Like Foster II

        raise NotImplementedError('TODO')

    def seriesparallelRL(self):
        """(R('R1') | L('L1')) + (R('R2') | L('L2'))"""
        # Like Foster I

        raise NotImplementedError('TODO')

    def fosterI(self, lexpr):

        expr = lexpr.partfrac(combine_conjugates=True)

        net = None
        for term in expr.as_ordered_terms():
            net = series(net, self.parallelRLC(LaplaceDomainExpression(term)))
        return net

    def fosterII(self, lexpr):

        expr = (1 / lexpr).partfrac(combine_conjugates=True)

        net = None
        for term in expr.as_ordered_terms():
            net = parallel(net, self.seriesRLC(
                LaplaceDomainExpression(1 / term)))
        return net

    def cauerI(self, lexpr):

        coeffs = lexpr.continued_fraction_coeffs()

        net = None
        for m, coeff in enumerate(reversed(coeffs)):
            n = len(coeffs) - m - 1
            if n & 1 == 0:
                net = series(net, self.seriesRL(
                    LaplaceDomainExpression(coeff)))
            else:
                net = parallel(net, self.parallelGC(
                    LaplaceDomainExpression(1 / coeff)))
        return net

    def cauerII(self, lexpr):

        coeffs = (1 / lexpr).continued_fraction_inverse_coeffs()

        net = None
        for m, coeff in enumerate(reversed(coeffs)):
            n = len(coeffs) - m - 1
            if n & 1 == 0:
                if not (n == 0 and coeff == 0):
                    net = parallel(net, self.seriesRL(
                        LaplaceDomainExpression(1 / coeff)))
            else:
                net = series(net, self.parallelGC(
                    LaplaceDomainExpression(coeff)))
        return net

    def network(self, lexpr, form='default'):

        if form == 'default':
            form = 'cauerI'

        lexpr = expr(lexpr)
        if not lexpr.is_impedance:
            raise ValueError('Expression needs to be an impedance')
        lexpr = lexpr.laplace()

        # Should test if a positive real function.

        try:
            method = getattr(self, form)
        except AttributeError:
            raise ValueError(
                '''Unknown form %s, known forms include: cauerI, cauerII, fosterI, fosterII''' % form)

        net = method(lexpr)
        return net


def network(lexpr, form='default'):

    return Synthesis().network(lexpr, form)
