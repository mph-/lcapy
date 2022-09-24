"""This module provides the Assumptions class.  It manages
assumptions for expressions needed for inverse Laplace transforms
and expression simplification.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from copy import copy
from .sym import ssym
from .acdc import is_dc, is_ac, is_causal


class Assumptions(dict):
    """Manage expression assumptions.  These include 'ac', 'dc', 'causal'
    and 'unknown' to specify expression behavior for t < 0, 'positive',
    'complex', etc., to help with expression simplification, 'nid' as a noise
    identifier, and `omega` for the phasor angular frequency.

    'ac', 'dc', 'causal' and 'unknown' refer to the time domain
    expression.  They can be inferred for a time domain expression but
    can be overridden.  The inference is performed lazily.  In
    particular, it is performed when converting from the time domain
    to the Laplace or Phasor domains.  The latter is required to
    handle the case of conversions from the Phasor to the Laplace
    domain.  Conversions from the Fourier domain to the Laplace domain
    go via the time domain and so the inference is performed while in
    the time domain.

    """

    def __init__(self, *args, **kwargs):

        super(Assumptions, self).__init__(*args)

        for assumption, value in kwargs.items():
            self.set(assumption, value)

    def set(self, assumption, value):

        if assumption in ('dc', 'ac', 'causal', 'unknown'):
            if value:
                self.pop('dc', None)
                self.pop('ac', None)
                self.pop('causal', None)
                self.pop('unknown', None)
                self[assumption] = value
            else:
                self.pop(assumption, None)
        else:
            self[assumption] = value

    def get(self, assumption, default=False):

        return super(Assumptions, self).get(assumption, default)

    def merge(self, **assumptions):

        new = self.copy()
        # TODO: warn if have mutually exclusive assumptions.
        # Note, the last one overrides.
        for assumption, value in assumptions.items():
            new.set(assumption, value)
        return new

    def copy(self):

        return copy(self)

    def sympy_assumptions(self):
        """Return dict of the SymPy assumptions such as complex, positive, etc."""

        assumptions = {}
        for assumption, value in self.items():
            if assumption not in ('nid', 'ac', 'dc', 'causal', 'unknown'):
                assumptions[assumption] = value
        return assumptions

    def infer_from_expr(self, expr):

        if expr.is_transform_domain:
            self.set('unknown', True)
            return

        var = expr.var
        if is_dc(expr, var):
            self.set('dc', True)
            return

        if is_ac(expr, var):
            self.set('ac', True)
            return

        if is_causal(expr, var):
            self.set('causal', True)
            return

        self.set('unknown', True)

    @property
    def has_unspecified(self):

        if 'ac' in self:
            return False
        if 'dc' in self:
            return False
        if 'causal' in self:
            return False
        if 'unknown' in self:
            return False
        return True

    def merge_and_infer(self, expr, **assumptions):
        """Override assumptions with specified assumptions.
        If none specified, infer them."""

        assumptions = self.merge(**assumptions)
        if assumptions.has_unspecified:
            assumptions.infer_from_expr(expr)
        return assumptions

    @property
    def is_ac(self):

        return self.get('ac')

    @property
    def is_dc(self):

        return self.get('dc')

    @property
    def is_causal(self):

        return self.get('causal')

    @property
    def is_unknown(self):

        return self.get('unknown')

    def convolve(self, expr, x):

        # The assumptions refer to the time-domain signal or
        # impulse-response.  Thus we want to propagate the assumptions
        # for the voltage or current signal.  The causal assumptions
        # are only required for s-domain expressions.
        #
        # The ac and dc assumptions were used for inverse Laplace
        # transforms but this sometimes gave bogus results.  For
        # example, when applying the eigenfunction for convolution,
        # `exp(q * t)`, to an LTI system, an additional transient
        # response is produced.  The transient response was not
        # removed and thus the result was not ac.

        # A similar error occured when DC was applied to an LTI system
        # since the inverse Laplace transform also produces a
        # transient response.  For example, let `x(t) = c` and `H(s) =
        # a / (s + a)`.  In this example, `X(s) = c / s$` and `ILT(X *
        # H)` yields `c - c * exp(-a * t)` for `t >= 0`.  Here the
        # second term is an artefact of the unilateral Laplace
        # transform since `x(t)` is treated as `x(t) u(t)`.  The true
        # result using time-domain convolution assuming `H(s)` is
        # causal is `c`.

        # The original signal can only be recovered from an inverse
        # unilateral Laplace transform if we have the region of
        # convergence.

        assumptions = self.copy()

        # If x is constant then copy assumptions.
        if x.var is None:
            return assumptions

        assumptions.set('unknown', True)

        if self.is_causal:
            if x.is_causal:
                assumptions.set('causal', True)
            elif expr.var == ssym:
                expr2 = x.sympy
                # Assume H(s) * s**n is causal if H(s) is causal.
                if expr2 == ssym or (expr2.is_Pow and expr2.args[0] == ssym):
                    assumptions.set('causal', True)
        return assumptions

    def add(self, x):

        assumptions = self.copy()
        assumptions.set('unknown', True)

        if self.is_causal and x.is_causal:
            assumptions.set('causal', True)
        elif self.is_dc and x.is_dc:
            assumptions.set('dc', True)
        elif False and self.is_ac and x.is_ac:
            # Hmmm, need to check if same frequency.
            assumptions.set('ac', True)
        return assumptions
