"""This module handles sequences.

Copyright 2020 Michael Hayes, UCECE

"""

from .expr import ExprList, expr

# Perhaps subclass numpy ndarray?  But then could not have symbolic
# elements in the sequence.

class Sequence(ExprList):

    def __init__(self, seq, n=None, evaluate=False, var=None):

        super (Sequence, self).__init__(seq, evaluate)

        if n is None:
            n = list(range(len(seq)))
        
        self.n = n
        self.var = var

    @property
    def vals(self):
        return list(self)
        
    def __getitem__(self, n):
        """Note this returns the element with index matching n.
        This is not necessarily the nth element in the sequence."""

        # TODO, support slices, etc.
        try:
            nindex = self.n.index(n)
        except ValueError:
            return expr(0)            

        return super(Sequence, self).__getitem__(nindex)

    def prune(self):
        """Remove zeros from ends of sequence.

        {0, 0, 1, 2, 3, 0} -> {1, 2, 3}"""

        vals = self.vals
        
        m1 = 0        
        while vals[m1] == 0:
            m1 += 1

        m2 = len(vals) - 1
        while vals[m2] != 0:
            return Sequence(vals[m1:], self.n[m1:])
        
        while vals[m2] == 0:
            m2 -= 1        
        return Sequence(vals[m1:m2 + 1], self.n[m1:m2 + 1])

    def zeropad(self, M):
        """Add M zeros to end of sequence:
        
        For example, with M = 3

        {1, 2, 3} -> {1, 2, 3, 0, 0, 0}"""

        vals = self.vals
        n = self.n

        zero = expr(0)
        for m in range(M):
            vals.append(zero)

        n = self.n + list(range(self.n[-1] + 1, len(vals)))
        return self.__class__(vals, n=n, var=self.var)        
    
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

    def plot(self, n=None, **kwargs):
        """Plot the sequence.  If n is not specified, it defaults to the
        range (-20, 20).  n can be a vector of specified instants, a
        tuple specifing the range, or a constant specifying the
        maximum value with the minimum value set to 0.

        kwargs include:
        axes - the plot axes to use otherwise a new figure is created
        xlabel - the x-axis label
        ylabel - the y-axis label
        xscale - the x-axis scaling, say for plotting as ms
        yscale - the y-axis scaling, say for plotting mV
        in addition to those supported by the matplotlib plot command.
        
        The plot axes are returned."""

        return self.as_impulses().plot(n, **kwargs)

    def __call__(self, arg, **assumptions):
        """Transform domain or substitute arg for variable. 
        
        Substitution is performed if arg is a tuple, list, numpy
        array, or constant.  If arg is a tuple or list return a list.
        If arg is an numpy array, return numpy array.

        Domain transformation is performed if arg is a domain variable
        or an expression of a domain variable.

        See also evaluate.

        """

        return self.as_impulses().__call__(arg, **assumptions)

    def ZT(self, **assumptions):
        return self.as_impulses().ZT(**assumptions)    

    def _repr_pretty_(self, p, cycle):
        """This is used by jupyter notebooks to display an expression using
        unicode.  It is also called by IPython when displaying an
        expression.""" 

        # Note, the method in ExprPrint is bypassed since list
        # has this methodx.
        
        from .printing import pretty
        p.text(self.pretty())

    def lfilter(self, b=[], a=[1]):

        x = self.vals
        y = []

        a0 = a[0]

        for n, x1 in enumerate(x):
            y.append(expr(0))
            
            for m, b1 in enumerate(b):
                try:
                    y[-1] += b1 * x[n - m] / a0
                except:
                    pass                

            yn = y[-1]
            for m, a1 in enumerate(a[1:]):
                try:
                    yn += a1 * y[-m - 2] / a0
                except:
                    pass
            y[-1] = yn
                
        #n = self.n + list(range(self.n[-1] + 1, len(y)))
        n = self.n
        return self.__class__(y, n=n, var=self.var)
    
