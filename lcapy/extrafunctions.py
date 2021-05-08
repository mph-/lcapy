import sympy as sym
from sympy.core import S, Integer
from sympy.core.logic import fuzzy_not

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
    def eval(cls, nval):
        """
        Evaluates the discrete unit step function.
        """
        
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

        
class sincn(sym.Function):
    
    @classmethod
    def eval(cls, val):
        """
        Evaluates the normalised sinc (cardinal sine) function.
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
        Evaluates the unnormalised sinc (cardinal sine) function.
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

