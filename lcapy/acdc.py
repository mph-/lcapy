from lcapy.sympify import sympify1
import sympy as sym
from sympy import symbols, I, exp, cos, pi, sin

class Causal(object):

    def _has_causal_factor(self, expr):

        factors = expr.as_ordered_factors()
        for factor in factors:
            if factor == 0:
                return True
            if (not factor.is_Function 
                or factor.func not in (sym.Heaviside, sym.DiracDelta)):
                continue

            p = sym.Poly(factor.args[0], self.var)
            coeffs = p.all_coeffs()
            if len(coeffs) != 2:
                return False

            if (coeffs[0].is_positive 
                and (coeffs[1].is_positive or coeffs[1].is_zero)):
                return True

        return False        

    def _is_causal(self, expr):

        terms = expr.as_ordered_terms()
        for term in terms:
            if not self._has_causal_factor(term):
                return False
                
        return True

    def __init__(self, expr, var):        

        self.var = getattr(var, 'expr', sympify1(var)) 
        self.causal = self._is_causal(getattr(expr, 'expr', sympify1(expr)))


class DC(object):

    def _is_dc(self, expr):

        if self.var in expr.free_symbols:
            return False

        terms = expr.as_ordered_terms()
        for term in terms:
            for factor in term.as_ordered_factors():
                n, d = factor.as_numer_denom()
                if not ((n.is_Symbol or n.is_number) and (d.is_Symbol or d.is_number)):
                    return False
        return True

    def __init__(self, expr, var):        

        self.var = getattr(var, 'expr', sympify1(var))
        self.dc = self._is_dc(getattr(expr, 'expr', sympify1(expr)))


class AC(object):

    def _find_freq_phase(self, expr):

        self.freq = 0

        if expr.func == cos:
            self.phase = 0
        elif expr.func == sin:
            self.phase = pi / 2
        else:
            raise ValueError('%s not sin/cos' % expr)
            
        p = sym.Poly(expr.args[0], self.var)
        coeffs = p.all_coeffs()
        if len(coeffs) != 2:
            return False

        self.phase += coeffs[1]
        self.freq = coeffs[0] / (2 * pi)
        return True

    def _is_ac(self, expr):

        # Convert sum of exps into sin/cos
        expr = expr.rewrite(cos).combsimp().expand()
        
        factors = expr.as_ordered_factors()
            
        self.amp = 1
        for factor in factors:
            if factor.is_Function:
                if factor.func not in (cos, sin):
                    return False
                if not self._find_freq_phase(factor):
                    return False
            elif is_dc(factor, self.var):
                self.amp *= factor
            else:
                return False
        return True

    def __init__(self, expr, var):

        self.var = getattr(var, 'expr', sympify1(var))
        self.amp = 0
        self.freq = 0
        self.phase = 0
        self.ac = self._is_ac(getattr(expr, 'expr', sympify1(expr)))


def is_dc(expr, var):
    return DC(expr, var).dc


def is_ac(expr, var):
    return AC(expr, var).ac


def is_causal(expr, var):
    return Causal(expr, var).causal
