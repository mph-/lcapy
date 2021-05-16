"""This module provides discrete-time filter support.

Copyright 2021 Michael Hayes, UCECE

"""

from .expr import expr, equation
from .discretetime import n, z
import sympy as sym

def isiterable(arg):

    return hasattr(arg, '__iter__')


class DTFilter(object):

    def __init__(self, b, a):
        """Create discrete-time filter where `b` is a list or array of
        numerator coefficients and `a` is a list of array of
        denominator coefficients."""

        if not isiterable(b):
            b = (b, )
        if not isiterable(a):
            a = (a, )            
        
        self.a = a
        self.b = b

    def transfer_function(self):

        Nl = len(self.a)
        Nr = len(self.b)
        
        # numerator of H(z)
        num = 0 * z 
        for i in range(Nr): 
            num += self.b[i] * z**(-i)
            
        # denominator for H(z)
        denom = self.a[0] * z**0  
        for k in range(1, Nl):
            az = self.a[k] * z**(-k)
            denom += az
  
        # collect with respect to positive powers of the variable z
        num = sym.collect(sym.expand(num * z**Nl), z)
        denom = sym.collect(sym.expand(denom * z**Nl), z)
        
        Hz = expr(sym.simplify(num / denom))
        Hz.is_causal = True
        return Hz

    def impulse_response(self):

        H = self.transfer_function()
        return H(n)

    def difference_equation(self, input='x', output='y'):

        rhs = 0 * n

        for m, bn in enumerate(self.b):
            rhs += bn * expr('%s(n - %d)' % (input, m))

        for m, an in enumerate(self.a[1:]):
            rhs -= an * expr('%s(n - %d)' % (output, m + 1))

        lhs = expr('y(n)')
        e = equation(lhs, rhs)

        return e
    
