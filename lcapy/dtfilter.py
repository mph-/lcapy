"""This module provides discrete-time filter support.

Copyright 2021 Michael Hayes, UCECE

"""

from .expr import expr, equation
from .discretetime import n, z, seq
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
        """Return discrete-time transfer function."""                

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
        """Return discrete-time impulse response."""        

        H = self.transfer_function()
        return H(n)

    def difference_equation(self, input='x', output='y'):
        """Return difference equation."""

        rhs = 0 * n

        for m, bn in enumerate(self.b):
            rhs += bn * expr('%s(n - %d)' % (input, m))

        for m, an in enumerate(self.a[1:]):
            rhs -= an * expr('%s(n - %d)' % (output, m + 1))

        lhs = expr('y(n)')
        e = equation(lhs, rhs)

        return e

    def zdomain_initial_response(self, ic):
        """Return zdomain response due to initial conditions."""        

        Nl = len(self.a)
        
        # Denominator for Yi(z)
        denom = self.a[0] * z**0  
        num = 0 * z
        for k in range(1, Nl):
            az = self.a[k] * z**(-k)
            denom += az
            # Numerator for Yi(z)
            y0 = 0 * z
            for i in range(0, k):
                y0 += ic[i] * z**(i + 1)
                num += az * y0
  
        # Collect with respect to positive powers of the variable z
        num = sym.collect(sym.expand(num * z**Nl), z)
        denom = sym.collect(sym.expand(denom * z**Nl), z)
  
        Yzi = expr(-sym.simplify(num / denom))
        Yzi.is_causal = True
        
        return Yzi

    def initial_response(self, ic):
        """Return response due to initial conditions."""

        Yzi = self.zdomain_initial_response(ic)
        return Yzi(n)
    
    def response(self, x, ic, ns):

        NO = len(ic)
        
        if NO !=  (len(self.a)-1):
            raise ValueError("Expected %d initial conditions, got %d" % (len(self.a) - 1, NO))
    
        # Number of n vals, including last one
        alln = list(range(ns[0], ns[1] + 1))
        Nn = len(alln)
  
        # Order right hand side
        Nr = len(self.b)
  
        y_tot = list(ic[-1::-1]) + Nn * [0]
  
        self.a_r = self.a[-1:-1-NO:-1]
        for i, nval in enumerate(alln):
            # Get previous y vals (sliding window)
            pre_y = y_tot[i:i + NO]
    
        # Calculate rhs of new value
        rhs = sum(self.b[l] * inp_func(nval - l) for l in range(Nr))
    
        # Add lhs
        y_tot[i + NO] = -1 / self.a[0] * sum(csi * ysi for csi, ysi in zip(self.a_r, pre_y)) + rhs
    
        # insert underscore at n = 0
        # the following does not work, if the result will be plotted with xn.plot()
        #try:
        #i = alln.index(0)
        #y_tot[i + NO] = '_(%s)'%(y_tot[i + NO])
        #except:
        #pass
  
        # solution, without initial values
        ret_seq = seq(y_tot[NO:])  
  
        return ret_seq
