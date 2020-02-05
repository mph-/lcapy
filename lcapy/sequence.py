"""This module handles sequences.

Copyright 2020 Michael Hayes, UCECE

"""

from .expr import ExprList

class Sequence(ExprList):

    def __init__(self, seq, n=None):

        super (Sequence, self).__init__(seq)

        # Save the indexes.  Ideally, should annotate which item
        # in sequence corresponds to n = 0.
        self.n = n

    def latex(self):

        items = []
        for v1, n1 in zip(self.n, self):
            s = v.latex()
            if n1 == 0:
                s = r'\underline{%s}' % v1
            items.append(s)

        return '\left{%s\right\}' % ', '.join(items)
    
            
