"""This module handles sequences.

Copyright 2020 Michael Hayes, UCECE

"""

from .expr import ExprList, expr

# Perhaps subclass numpy ndarray?  But then could not have symbolic
# elements in the sequence.  Perhaps have flavours for n-domain and
# k-domain sequences?

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

    def __str__(self):

        items = []
        for v1, n1 in zip(self, self.n):
            s = str(v1)
            
            if n1 == 0:
                s = '_' + s
            items.append(s)

        return r'{%s}' % ', '.join(items)
    
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
        for n1 in range(0, self.n[0]):
            if n1 == 0:
                s = r'\underline{0}'
            else:
                s = '0'
            items.append(s)                
            
        for v1, n1 in zip(self, self.n):
            try:
                s = v1.latex()
            except:
                s = str(v1)
            
            if n1 == 0:
                s = r'\underline{%s}' % v1
            items.append(s)

        return r'\left\{%s\right\}' % ', '.join(items)
    
    def pretty(self, **kwargs):

        # FIXME, this does not work well if values require multiple
        # lines, such as a^2
        # See _print_seq in pretty.py of SymPy but how to denote
        # the underline?
        
        items = []
        for n1 in range(0, self.n[0]):
            if n1 == 0:
                s = r'_0'
            else:
                s = '0'
            items.append(s)                            
        
        for v1, n1 in zip(self, self.n):
            try:
                s = v1.pretty()
            except:
                s = str(v1)
            
            if n1 == 0:
                s = '_' + s
            items.append(s)

        return r'{%s}' % ', '.join(items)
    
    def as_impulses(self, var=None):
        """Convert to discrete-time signal in the form of
        a weighted sum of delayed impulses.  For example,
        {1, 2, 3} -> ui[n] + 2 * ui[n - 1] + 3 * ui[n - 2]"""

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

        For example, Sequence([1, 1]).extent() = 2
                     Sequence([1, 0, 1]).extent() = 3
                     Sequence([0, 1, 0, 1]).extent() = 3
        """
        
        from numpy import argwhere

        # Note, each element is an Expr.
        nz = [elt != 0 for elt in self]
        w = argwhere(nz)
        if len(w) == 0:
            return 0

        return (max(w) - min(w))[0] + 1

    def discrete_time_fourier_transform(self, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""

        return self.as_impulses().discrete_time_fourier_transform(**assumptions)

    def DTFT(self, **assumptions):
        """Convert to Fourier domain using discrete time Fourier transform."""
    
        return self.discrete_time_fourier_transform(**assumptions)

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
        """Convert sequence to n-domain or k-domain expression.
        For example, seq((1, 2, 3))(n)"""

        from .nexpr import n
        from .kexpr import k

        if id(arg) in (id(n), id(k)):
            return self.as_impulses(arg)
        if arg in (n, k):
            return self.as_impulses(arg)        

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

    def lfilter(self, b=None, a=None):
        """Implement digital filter specified by a transfer function.  The
        transfer function is described by a vector `b` of coefficients
        for the numerator and a `a` vector of coefficients for the
        denominator.

        For a FIR filter a = [1]."""

        if b is None:
            b = []
        if a is None:
            a = [1]
        
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
    
    def convolve(self, h, mode='full'):
        """Convolve with h."""

        x = self
        h = Sequence(h)
        
        Lx = x.extent()
        Lh = h.extent()
        Ly = Lx + Lh - 1
        
        if mode == 'full':
            x = x.zeropad(Ly - Lx)
        elif mode == 'same':
            x = x.zeropad(max(Lx, Ly) - Lx)            
        else:
            raise ValueError('Unknown mode ' + mode)
        
        return x.lfilter(h, a=[1])
    
