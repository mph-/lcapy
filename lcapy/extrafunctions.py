"""This module defines additional functions to those defined by SymPy.
They are for internal use by Lcapy.  The user versions are defined in
function.py.

Copyright 2020--2021 Michael Hayes, UCECE

"""

import sympy as sym
from sympy.core import S, Integer
from sympy.core.logic import fuzzy_not
from .config import unitstep_zero

class UnitImpulse(sym.Function):

    is_integer = True
    
    @classmethod
    def eval(cls, nval):
        """
        Evaluates the discrete unit impulse function.
        """
        
        if nval.is_zero:
            return S.One
        elif fuzzy_not(nval.is_zero):
            return S.Zero


class UnitStep(sym.Function):

    is_integer = True
    
    @classmethod
    def eval(cls, nval, zero=None):
        """
        Evaluates the discrete unit step function.   This is defined
        as 1 for n >= 0 and 0 otherwise.  
        """

        if nval.is_zero:
            if zero is None:
                zero = unitstep_zero
            return zero
        
        if nval.is_nonnegative:
            return S.One
        elif nval.is_negative:
            return S.Zero

        
class rect(sym.Function):

    @classmethod
    def eval(cls, val):
        """
        Evaluates the rectangle function.
        """

        if val.is_Number:
            if val < -0.5 or val > 0.5:
                return S.Zero
            return S.One

    def rewrite(self, *args, **hints):

        x = self.args[0]
        return sym.Heaviside(x + S.Half) - sym.Heaviside(x - S.Half)


class dtrect(sym.Function):

    @classmethod
    def eval(cls, val):
        """
        Evaluates the rectangle function for discrete-time signals.
        """

        if val.is_Number:
            if val < -0.5 or val >= 0.5:
                return S.Zero
            return S.One

    def rewrite(self, *args, **hints):

        x = self.args[0]
        return UnitStep(x + S.Half) - UnitStep(x - S.Half)

    
class dtsign(sym.Function):

    @classmethod
    def eval(cls, val):
        """
        Evaluates the signum function for discrete-time signals.
        """

        if val.is_Number:
            if val >= 0:
                return S.One
            return S.Zero

        
class sincn(sym.Function):
    
    @classmethod
    def eval(cls, val):
        """
        Evaluates the normalized sinc (cardinal sine) function.
        This is what NumPy uses but not SymPy.
        """

        if val.is_Number:
            if val == 0:
                return S.One
            x = sym.pi * val
            return sym.sin(x) / x
        
    def rewrite(self, *args, **hints):

        x = sym.pi * self.args[0]
        return sym.sin(x) / x


class sincu(sym.Function):
    
    @classmethod
    def eval(cls, val):
        """
        Evaluates the unnormalized sinc (cardinal sine) function.
        This is what SymPy uses but not NumPy.
        """

        if val.is_Number:
            if val == 0:
                return S.One
            return sym.sin(val) / val
        
    def rewrite(self, *args, **hints):

        x = self.args[0]
        return sym.sin(x) / x    

class psinc(sym.Function):
    
    @classmethod
    def eval(cls, M, val):
        """
        Evaluates the periodic sinc function.
        """
        if val.is_Number:
            if val == 0:
                return S.One
            x = sym.pi * val
            return sym.sin(M * x) / (M * sym.sin(x))
        
    def rewrite(self, *args, **hints):

        M = self.args[0]        
        x = sym.pi * self.args[1]
        return sym.sin(M * x) / (M * sym.sin(x))

    
class tri(sym.Function):

    @classmethod
    def eval(cls, val):
        """
        Evaluates the triangle function.
        """

        if val.is_Number:
            
            if val >= 1:
                return S.Zero
            elif val <= -1:
                return S.Zero
            else:
                return 1 - abs(val)


class trap(sym.Function):

    @classmethod
    def eval(cls, val, alpha):
        """
        Evaluates the trapezoid function.   This is rect(t / alpha).convolve(rect(t)).
        trap(t, 0) = rect(t) and trap(t, 1) = tri(t).
        """

        if val.is_Number and alpha.is_number:
            absval = abs(val)
            foo = absval - 0.5

            if alpha == S.Zero:
                if val < -0.5 or val > 0.5:
                    return S.Zero
                return S.One
            
            if foo >= 0.5 * alpha:
                return S.Zero
            elif foo <= -0.5 * alpha:
                return S.One
            else:
                return 0.5 - foo / alpha

