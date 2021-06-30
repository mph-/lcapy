"""This module provides the SequenceExpression class to provide
common methods for the discrete-time and discrete-frequency
expressions.

Copyright 2020--2021 Michael Hayes, UCECE

"""

from .dexpr import DiscreteExpression
from .sequence import Sequence
from .functions import Heaviside, DiracDelta, rect, sign
from .extrafunctions import UnitStep, UnitImpulse, dtrect, dtsign
from numpy import arange


class SequenceExpression(DiscreteExpression):
    """Superclass of discrete-time and discrete-frequency expressions."""

    def __init__(self, val, **assumptions):

        super(SequenceExpression, self).__init__(val, **assumptions)

        mapping = {Heaviside: UnitStep,
                   DiracDelta: UnitImpulse,
                   rect: dtrect,
                   sign: dtsign}

        # Use discrete-time function variants, see also functions.py
        for old, new in mapping.items():
            if self.has(old):
                self.expr = self.replace(old, new).expr                
                
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

        start_trunc = False
        end_trunc = False
        
        if ni is None:
            n1 = self.first_index(ni)
            n2 = self.last_index(ni)        
            ni = arange(n1, n2 + 1)
            # Should search to handle cases such as 1, 1, 0, 1 but
            # when to stop?
            if self(n1 - 1) != 0:
                print('Warning: sequence truncated at n1=%d' % n1)
                start_trunc = True
            if self(n2 + 1) != 0:
                print('Warning: sequence truncated at n2=%d' % n2)
                end_trunc = True
            
        elif isinstance(ni, tuple):
            ni = arange(ni[0], ni[-1] + 1)            
            
        v = self(ni)
        
        return self.seqcls(v, ni, evaluate=evaluate,
                           start_trunc=start_trunc, end_trunc=end_trunc)
       
