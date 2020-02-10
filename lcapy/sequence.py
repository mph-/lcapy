"""This module handles sequences.

Copyright 2020 Michael Hayes, UCECE

"""

from .expr import ExprList

# Perhaps subclass numpy ndarray?  But then could not have symbolic
# elements in the sequence.

class Sequence(ExprList):

    def __init__(self, seq, n=None, evaluate=False, var=None):

        super (Sequence, self).__init__(seq, evaluate)

        if n is None:
            n = list(range(len(seq)))
        
        self.n = n
        self.var = var

    def latex(self):

        items = []
        for v1, n1 in zip(self, self.n):
            try:
                s = v1.latex()
            except:
                s = str(v1)
            
            if n1 == 0:
                s = r'\underline{%s}' % v1
            items.append(s)

        return r'\left\{%s\right\}' % ', '.join(items)
    
    def pretty(self):

        items = []
        for v1, n1 in zip(self, self.n):
            try:
                s = v1.pretty()
            except:
                s = str(v1)
            
            if n1 == 0:
                s = '_%s' % v1
            items.append(s)

        return r'{%s}' % ', '.join(items)
    
    def as_impulses(self, var=None):

        from .functions import unitimpulse
        from sympy import Add

        if var is None:
            var = self.var
        if var is None:
            raise ValueError('var not specified')

        # This can reorder terms
        result = var * 0
        for v1, n1 in zip(self, self.n):        
            result += v1 * unitimpulse(var - n1)

        # but so does the unevaluated Add...
        # items = []
        # for v1, n1 in zip(self, self.n):        
        #     items.append(v1 * unitimpulse(var.var - n1))
        #
        # expr = Add(*items, evaluate=False)
        # result = var.__class__(expr)
        
        return result

    def extent(self):
        """Determine extent of the sequence.

        For example, Sequence([1, 1]).extent = 2
                     Sequence([1, 0, 1]).extent = 3
                     Sequence([0, 1, 0, 1]).extent = 3
        """
        
        from numpy import argwhere

        # Note, each element is an Expr.
        nz = [elt != 0 for elt in self]
        w = argwhere(nz)
        if len(w) == 0:
            return 0

        return (max(w) - min(w))[0] + 1
