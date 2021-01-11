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

        
class sinc(sym.Function):

    @classmethod
    def eval(cls, val):
        """
        Evaluates the sinc (cardinal sine) function.
        """

        if val.is_Number:
            if val == 0:
                return S.One
            x = sym.pi * val
            return sym.sin(x) / x
        
    def rewrite(self, *args, **hints):

        x = sym.pi * self.args[0]
        return sym.sin(x) / x
    
