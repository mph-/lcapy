"""This module handles sequences.

Copyright 2020 Michael Hayes, UCECE

"""

from .expr import ExprList


class Sequence(ExprList):

    def __init__(self, seq, n=None, evalf=False):

        super (Sequence, self).__init__(seq, evalf)
        self.n = n

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
                s = '_%s_' % v1
            items.append(s)

        return r'{%s}' % ', '.join(items)
    
    def sum_delayed_impulses(self, var):

        from .discretetime import unitimpulse
        
        result = var * 0
        for v1, n1 in zip(self, self.n):        
            result += v1 * unitimpulse(var - n1)
        return result


    
