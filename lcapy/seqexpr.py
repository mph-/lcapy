"""This module provides the SequenceExpression class to provide
common methods for the discrete-time and discrete-frequency
expressions.

Copyright 2020--2022 Michael Hayes, UCECE

"""

from .dexpr import DiscreteExpression
from .sequence import Sequence
from .functions import function_mapping
from .extrafunctions import UnitStep, UnitImpulse, dtrect, dtsign
from sympy import Heaviside
from warnings import warn


class SequenceExpression(DiscreteExpression):
    """Superclass of discrete-time and discrete-frequency expressions."""

    def __init__(self, val, **assumptions):

        super(SequenceExpression, self).__init__(val, **assumptions)

        def remap_continuous_discrete(expr):

            def query(expr):
                return expr.is_Function and expr.func in function_mapping

            def value(expr):
                # Note, Heaviside behaviour changed in SymPy-1.9
                # The second arg now defaults to 1/2, so we need
                # to remove it.
                if expr.func == Heaviside:
                    return UnitStep(expr.args[0])
                return function_mapping[expr.func](*expr.args)

            # Use discrete-time function variants, see also functions.py
            return expr.replace(query, value)

        self.expr = remap_continuous_discrete(self.expr)

    def first_index(self, ni=None):

        if ni is None:
            ni = (-10, 10)
        if isinstance(ni, tuple):
            ni = range(*ni)

        # Desire SymPy equivalent to argmax
        n2 = ni[-1]
        for n in reversed(ni):
            if self(n) != 0:
                n2 = n
        return n2

    def last_index(self, ni=None):

        if ni is None:
            ni = (-10, 10)
        if isinstance(ni, tuple):
            ni = range(*ni)

        # Desire SymPy equivalent to argmin
        n1 = ni[0]
        for n in ni:
            if self(n) != 0:
                n1 = n
        return n1

    def seq(self, ni=None, evaluate=False):
        """Create a Sequence.

        >>> a = x.seq()

        The sequence indices are specified with the optional `ni` argument.
        For example:

        >>> a = x.seq(ni=(-1, 0, 1, 2))

        If the `ni` argument is not specified, the sequence indices
        are enumerated from 0.

        The sequence indices can be found using the `n` attribute.
        This returns a list.

        >>> a = x.seq().n
        [-1, 0, 1, 2]
        """

        from numpy import arange

        start_trunc = False
        end_trunc = False

        if ni is None:
            n1 = self.first_index(ni)
            n2 = self.last_index(ni)
            ni = arange(n1, n2 + 1)
            # Should search to handle cases such as 1, 1, 0, 1 but
            # when to stop?
            if self(n1 - 1) != 0:
                warn('Sequence truncated at n1=%d' % n1)
                start_trunc = True
            if self(n2 + 1) != 0:
                warn('Sequence truncated at n2=%d' % n2)
                end_trunc = True

        elif isinstance(ni, tuple):
            ni = arange(ni[0], ni[-1] + 1)

        v = self(ni)

        return self.seqcls(v, ni, evaluate=evaluate,
                           start_trunc=start_trunc, end_trunc=end_trunc)

    @property
    def is_stable(self):
        """Return True if all the poles of the signal's z-transform
        are within the unit circle.

        See also is_marginally_stable."""

        poles = self.poles(aslist=True)
        for pole in poles:
            # If poles not constant, assume they are inside the unit circle.
            if pole.is_constant and abs(pole).fval >= 1:
                return False
        return True

    @property
    def is_marginally_stable(self):
        """Return True if all the poles of the signal's z-transform
        are within or on the unit circle.

        See also is_stable."""

        poles = self.poles(aslist=True)
        for pole in poles:
            # If poles not constant, assume they are inside the unit circle.
            if pole.is_constant and abs(pole).fval > 1:
                return False
        return True
