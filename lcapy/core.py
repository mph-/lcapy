"""
This module provides the core functions and classes for Lcapy.

To print the rational functions in canonical form (with the highest
power of s in the denominator with a unity coefficient), use
print(x.canonical()) or x.pprint() for pretty printing.

For additional documentation, see the Lcapy tutorial.

Copyright 2014 Michael Hayes, UCECE
"""

from __future__ import division
from warnings import warn
import numpy as np
import sympy as sym
from sympy.utilities.lambdify import lambdify
# Note imports at bottom to avoid circuit dependencies


__all__ = ('pprint', 'pretty', 'latex', 'DeltaWye', 'WyeDelta', 'tf', 
           'zp2tf', 'poles', 'zeros', 'residue', 'residues', 'partfrac',
           'general', 'canonical', 'ZPK', 'inverse_laplace', 'initial_value',
           'transient_response', 'response', 'final_value', 's', 'sExpr', 
           't', 'tExpr', 'cExpr')

class Expr(object):

    @property
    def expr(self):    
        return self.val
    
    
    def __init__(self, val, simplify=True, real=False):
        
        if isinstance(val, sExpr):
            val = val.val
            
        if real and isinstance(val, str) and val.isalnum() and not val[0].isdigit():
            val = sym.symbols(val, real=True)
        val = sym.sympify(val)

        if simplify:
            val = val.cancel()

        self.val = val


    def __str__(self):

        return self.val.__str__()


    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.val)

    
    def __abs__(self):
        """Absolute value"""
        
        return self.__class__(abs(self.val))


    def __neg__(self):
        """Negation"""
        
        return self.__class__(-self.val)


    def __rdiv__(self, x):
        """Reverse divide"""

        x = self.__class__(x)
        return self.__class__(x.val / self.val)


    def __rtruediv__(self, x):
        """Reverse true divide"""
            
        x = self.__class__(x)
        return self.__class__(x.val / self.val)


    def __mul__(self, x):
        """Multiply"""
        
        x = self.__class__(x)
        return self.__class__(self.val * x.val)


    def __rmul__(self, x):
        """Reverse multiply"""
        
        x = self.__class__(x)
        return self.__class__(self.val * x.val)


    def __div__(self, x):
        """Divide"""

        x = self.__class__(x)
        return self.__class__(self.val / x.val)


    def __truediv__(self, x):
        """True divide"""

        x = self.__class__(x)
        return self.__class__(self.val / x.val)
    

    def __add__(self, x):
        """Add"""
        
        x = self.__class__(x)
        return self.__class__(self.val + x.val)
    

    def __radd__(self, x):
        """Reverse add"""
        
        x = self.__class__(x)
        return self.__class__(self.val + x.val)
    
    
    def __rsub__(self, x):
        """Reverse subtract"""
        
        x = self.__class__(x)
        return self.__class__(x.val - self.val)
    

    def __sub__(self, x):
        """Subtract"""
        
        x = self.__class__(x)
        return self.__class__(self.val - x.val)
    
    
    def __pow__(self, x):
        """Pow"""
        
        x = self.__class__(x)
        return self.__class__(self.val ** x.val)


    def __or__(self, x):
        """Parallel combination"""
        
        return self.parallel(x)
    
    
    def __eq__(self, x):
        """Equality"""

        if x == None:
            return False

        x = self.__class__(x)
        return self.val == x.val


    def __ne__(self, x):
        """Inequality"""

        if x == None:
            return True

        x = self.__class__(x)
        return self.val != x.val


    def parallel(self, x):
        """Parallel combination"""
        
        x = self.__class__(x)
        return self.__class__(self.val * x.val / (self.val + x.val))
    

    def pprint(self):
        """Pretty print"""
        print(self.pretty())


    def pretty(self):
        """Make pretty string"""
        return sym.pretty(self.val)


    def _pretty(self, arg):
        """Make pretty string"""
        return sym.pretty(self.val)


    def latex(self):
        """Make latex string"""

        return sym.latex(self.val)

    
    @property
    def is_number(self):

        return self.expr.is_number


    @property
    def is_constant(self):

        return self.expr.is_constant()



class sExpr(Expr):
    """s-domain expression or symbol"""
    
    var = sym.symbols('s')

    def differentiate(self):
        """Differentiate (multiply by s)"""
        
        return self.__class__(self.val * self.var)
    

    def integrate(self):
        """Integrate (divide by s)"""
        
        return self.__class__(self.val / self.var)
    
    
    def delay(self, T):
        """Apply delay of T seconds by multiplying by exp(-s T)"""
        
        T = self.__class__(T)
        return self.__class__(self.val * sym.exp(-s * T))


    def subs(self, *args, **kwargs):
        """Substitute symbol for s in expression"""

        return self.expr.subs(s, *args, **kwargs)        


    def omega(self):
        """Return expression with s = j omega"""

        omega = sym.symbols('omega')
        return self.subs(sym.I * omega)        


    def zeros(self):
        """Return zeroes of expression"""
        
        return zeros(self.expr, self.var)


    def poles(self):
        """Return poles of expression"""
        
        return poles(self.expr, self.var)


    def residues(self):
        """Return residues of expression"""
        
        return residues(self.expr, self.var)


    def canonical(self):
        """Convert rational function to canonical form with unity
        highest power of denominator.

        See also general, partfrac, and ZPK"""
        
        return self.__class__(canonical(self.expr, self.var), simplify=False)


    def general(self):
        """Convert rational function to general form

        See also canonical, partfrac, and ZPK"""
        
        return self.__class__(general(self.expr, self.var), simplify=False)


    def partfrac(self):
        """Convert rational function into partial fraction form.

        See also canonical, general, and ZPK"""

        return self.__class__(partfrac(self.expr, self.var), simplify=False)


    def ZPK(self):
        """Convert to pole-zero-gain (PZK) form.
        
        See also canonical, general, and partfrac"""
        
        return self.__class__(ZPK(self.expr, self.var), simplify=False)


    def initial_value(self):
        """Determine value at t = 0"""
        
        return initial_value(self.expr, self.var)


    def final_value(self):
        """Determine value at t = oo"""
        
        return final_value(self.expr, self.var)


    def _as_ratfun_parts(self):
        
        return _as_ratfun_parts(self.expr, self.var)


    def split_strictly_proper(self):
        
        N, D = self._as_ratfun_parts()

        Q, M = N.div(D)

        return Q.as_expr(), M / D


    def simplify(self):
        """Simplify"""

        val = sym.simplify(self.val)
        return self.__class__(val)
    
    
    def inverse_laplace(self):
        """Attempt inverse Laplace transform"""
        
        return inverse_laplace(self.expr, t, self.var)


    def transient_response(self, t=None):
        """Evaluate transient (impulse) response"""

        return transient_response(self.expr, t, self.var)
    
    
    def impulse_response(self, t=None):
        """Evaluate transient (impulse) response"""
        
        return self.transient_response(t)


    def step_response(self, t=None):
        """Evaluate step response"""
        
        return (self / self.var).transient_response(t)
    
    
    def frequency_response(self, f=None):
        """Evaluate frequency response"""
        
        fsym = sym.symbols('f')

        expr = self.val.subs(s, sym.I * 2 * sym.pi * fsym)
        
        if f == None:
            return expr
        
        func = lambdify(fsym, expr, modules="numpy")
        return np.array([func(f1) for f1 in f])


    def response(self, x, t):
        """Evaluate response to input signal x"""

        return response(self.expr, x, t, self.var)


    def decompose(self):

        ratfun, delay = _as_ratfun_delay(self.expr, self.var)
        
        N, D = _as_ratfun_parts(ratfun, self.var)

        return N, D, delay


class tExpr(Expr):
    """t-domain expression or symbol"""

    var = sym.symbols('t')


    def laplace(self):
        """Attempt Laplace transform"""
        
        F, a, cond = sym.laplace_transform(self.expr, self.var, s)
        return sExpr(F)


class cExpr(Expr):
    """Constant expression or symbol"""

    def __init__(self, val):
        
        # FIXME for expressions of cExpr
        if False and not isinstance(val, (cExpr, int, float, str)):
            raise ValueError('%s of type %s not int, float, or str' % (val, type(val)))

        super (cExpr, self).__init__(val, simplify=False, real=True)


s = sExpr('s')
t = tExpr('t')


def _guess_var(expr, var):

    if hasattr(expr, 'expr'):
        return expr.expr, expr.var

    if var != None:
        return expr, var

    #if not expr.is_rational_function():
    #    raise ValueError('Expression not a rational function')

    numer, denom = expr.as_numer_denom()
    try:
        P = sym.Poly(numer)
        if P.gens != ():
            return expr, P.gens[0]
    except:
        pass

    try:
        P = sym.Poly(denom)
        if P.gens != ():
            return expr, P.gens[0]
    except:
        raise ValueError('Cannot determine polynomial variable')


def pprint(expr):

    print(pretty(expr))


def pretty(expr):

    if hasattr(expr, 'pretty'):
        return expr.pretty()
    else:
        return sym.pretty(expr)


def latex(expr):

    if hasattr(expr, 'latex'):
        return expr.latex()
    else:
        return sym.latex(expr)


def DeltaWye(Z1, Z2, Z3):

    ZZ = (Z1 * Z2 + Z2 * Z3 + Z3 * Z1)
    return (ZZ / Z1, ZZ / Z2, ZZ / Z3)


def WyeDelta(Z1, Z2, Z3):

    ZZ = Z1 + Z2 + Z3
    return (Z2 * Z3 / ZZ, Z1 * Z3 / ZZ, Z1 * Z2 / ZZ)


def tf(numer, denom=1, var=None):
    """Create a transfer function from lists of the coefficient
    for the numerator and denominator"""

    if var == None:
        var = sym.symbols('s')

    N = sym.Poly(numer, var)
    D = sym.Poly(denom, var)

    return sExpr(N / D)


def _zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists of zeros and poles, 
    and from a constant gain"""
    
    if var == None:
        var = sym.symbols('s')

    zeros = sym.sympify(zeros)
    poles = sym.sympify(poles)

    if isinstance(zeros, (tuple, list)):
        zz = [(var - z) for z in zeros]
    else:
        zz = [(var - z) ** zeros[z] for z in zeros]

    if isinstance(zeros, (tuple, list)):
        pp = [1 / (var - p) for p in poles]
    else:
        pp = [1 / (var - p) ** poles[p] for p in poles]
        
    return sym.Mul(K, *(zz + pp))


def zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists of zeros and poles, 
    and from a constant gain"""

    return sExpr(_zp2tf(zeros, poles, K, var))


def poles(expr, var=None):

    expr, var = _guess_var(expr, var)

    numer, denom = expr.as_numer_denom()
    poles = sym.roots(sym.Poly(denom, var))
    return poles


def zeros(expr, var=None):

    expr, var = _guess_var(expr, var)

    numer, denom = expr.as_numer_denom()
    zeros = sym.roots(sym.Poly(numer, var))
    return zeros


def residue(expr, var, pole, poles):

    # Remove pole from list of poles; sym.cancel
    # doesn't always work, for example, for complex poles.
    poles2 = poles.copy()
    poles2[pole] -= 1

    numer, denom = expr.as_numer_denom()
    D = sym.Poly(denom, var)
    K = D.LC()

    D = [(var - p) ** poles2[p] for p in poles2]
    denom = sym.Mul(K, *D)

    d = sym.limit(denom, var, pole)

    if d != 0:
        tmp = numer / denom
        return sym.limit(tmp, var, pole)

    print("Trying l'hopital's rule")
    tmp = numer / denom
    tmp = sym.diff(tmp, var)
    
    return sym.limit(tmp, var, pole)


def residues(expr, var=None):

    expr, var = _guess_var(expr, var)

    N, D = _as_ratfun_parts(expr, var)

    Q, M = N.div(D)

    expr = M / D

    P = poles(expr, var)
    F = []
    R = []
    for p in P:
        
        # Number of occurrences of the pole.
        N = P[p]

        f = var - p

        if N == 1:
            F.append(f)
            R.append(residue(expr, var, p, P))
            continue

        # Handle repeated poles.
        expr2 = expr * f ** N
        for n in range(1, N + 1):
            m = N - n
            F.append(f ** n)
            R.append(sym.limit(sym.diff(expr2, var, m), var, p) / sym.factorial(m))

    return F, R, Q


def partfrac(expr, var=None):
    """Convert rational function into partial fraction form (standard form).

    See also canonical, general, and ZPK"""

    expr, var = _guess_var(expr, var)

    ratfun, delay = _as_ratfun_delay(expr, var)

    F, R, Q = residues(ratfun, var)

    expr = Q.as_expr()
    for f, r in zip(F, R):
        expr = expr + r / f

    if delay != 0:
        expr *= sym.exp(-var * delay)

    return expr


def _as_ratfun_parts(expr, var=None):
    
    expr, var = _guess_var(expr, var)

    if not expr.is_rational_function():
        raise ValueError('Expression not a rational function')

    numer, denom = expr.as_numer_denom()
    N = sym.Poly(numer, var)
    D = sym.Poly(denom, var)
    
    return N, D


def _as_ratfun_delay(expr, var=None):
    
    expr, var = _guess_var(expr, var)

    F = sym.factor(expr).as_ordered_factors()

    delay = sym.sympify(0)
    ratfun = sym.sympify(1)
    for f in F:
        b, e = f.as_base_exp()
        if b == sym.E and e.is_polynomial(var):
            p = sym.Poly(e, var)
            c = p.all_coeffs()
            if p.degree() == 1:
                delay -= c[0]
                if c[1] != 0:
                    ratfun *= sym.exp(c[1])
                continue

        ratfun *= f
            
    if not ratfun.is_rational_function():
        raise ValueError('Expression not a product of rational function and exponential')

    return ratfun, delay


def general(expr, var=None):
    """Convert rational function into general form.

    See also canonical, partfrac, and ZPK"""

    expr, var = _guess_var(expr, var)

    return sym.cancel(expr, var)



def canonical(expr, var=None):
    """Convert rational function into canonical form.

    See also general, partfrac, and ZPK"""

    expr, var = _guess_var(expr, var)

    ratfun, delay = _as_ratfun_delay(expr, var)

    N, D = _as_ratfun_parts(ratfun, var)

    K = sym.cancel(N.LC() / D.LC())
    if delay != 0:
        K = K * sym.exp(var * delay)

    # Divide by leading coefficient
    N = N.monic()
    D = D.monic()
    
    return K * (N / D)


def ZPK(expr, var=None):
    """Convert rational function into ZPK form.

    See also canonical, general, and partfrac"""

    expr, var = _guess_var(expr, var)

    ratfun, delay = _as_ratfun_delay(expr, var)

    N, D = _as_ratfun_parts(ratfun, var)

    K = sym.cancel(N.LC() / D.LC())
    if delay != 0:
        K = K * sym.exp(var * delay)

    zeros = sym.roots(N)
    poles = sym.roots(D)

    return _zp2tf(zeros, poles, K, var)


def _inverse_laplace(expr, t, var):

    ratfun, delay = _as_ratfun_delay(expr, var)

    N, D = _as_ratfun_parts(ratfun, var)

    Q, M = N.div(D)

    result1 = 0

    # Delayed time.
    td = t - delay

    if Q:
        C = Q.all_coeffs()
        for n, c in enumerate(C):
            result1 += c * sym.diff(sym.DiracDelta(td), t, len(C) - n - 1)
        
    expr = M / D

    P = poles(expr, var)
    result2 = 0
    for p in P:
        
        # Number of occurrences of the pole.
        N = P[p]

        if N == 0:
            continue

        f = var - p

        if N == 1:
            r = residue(expr, var, p, P)

            pc = p.conjugate()
            if pc != p and P.has_key(pc):
                # Remove conjugate from poles and process pole with its conjugate
                P[pc] = 0
                result2 += 2 * sym.re(r) * sym.exp(sym.re(p) * td) * sym.cos(sym.im(p) * td) + 2 * sym.im(r) * sym.exp(sym.re(p) * td) * sym.sin(sym.im(p) * td)
            else:
                result2 += r * sym.exp(p * td)
            continue

        # Handle repeated poles.
        expr2 = expr * f ** N
        for n in range(1, N + 1):
            m = N - n
            r = sym.limit(sym.diff(expr2, var, m), var, p) / sym.factorial(m)
            result2 += r * sym.exp(p * td) * td**(n - 1)

    return result1 + result2 * sym.Heaviside(td)


def inverse_laplace(expr, t=None, s=None):
    """Determine inverse Laplace transform of expression"""

    expr, s = _guess_var(expr, s)
    if t == None:
        t = sym.symbols('t')

    try:
        result = _inverse_laplace(expr, t, s)
    except:

        print('Determining inverse Laplace transform with sympy...')

        # Try splitting into partial fractions to help sympy.
        expr = partfrac(expr, s)

        # This barfs when needing to generate Dirac deltas
        from sympy.integrals.transforms import inverse_laplace_transform
        result = inverse_laplace_transform(expr, t, s)

    return tExpr(result)


def transient_response(expr, t=None, s=None):
    """Determine transient (impulse) response"""        

    if isinstance(t, sym.Expr):
        return inverse_laplace(expr, t, s)

    tv = sym.symbols('t')

    texpr = inverse_laplace(expr, tv, s)
    if t == None:
        return texpr

    print('Evaluating inverse Laplace transform...')
        
    func = lambdify(tv, texpr, ("numpy", "sympy", "math"))
        
    try:
        # FIXME for scalar t
        response = np.array([complex(func(t1)) for t1 in t])
        # The following does not work if all the sympy functions are not
        # converted to numpy functions.
        # response = func(t)
            
    except NameError:
        raise RuntimeError('Cannot evaluate inverse Laplace transform')

    except AttributeError:
        raise RuntimeError('Cannot evaluate inverse Laplace transform, probably have undefined symbols, such as Dirac delta')
        
    # The result should be real so quietly remove any imaginary component.
    return response.real


def response(expr, x, t=None, s=None):
    """Determine response to excitation x"""        

    expr, s = _guess_var(expr, s)

    ratfun, delay = _as_ratfun_delay(expr, s)
    N, D = _as_ratfun_parts(ratfun, s)
    Q, M = N.div(D)
    expr = M / D

    h = transient_response(expr, t, s)
    y = np.convolve(x, h)

    if Q:
        # Handle Dirac Deltas and derivatives.
        C = Q.all_coeffs()
        dt = np.diff(t)
        for n, c in enumerate(C):
            
            y += c * x

            x = np.diff(x) / dt
            x = np.hstack((x, 0))

    if delay != 0.0:
        td = t - delay

        import scipy.interpolate.interp1d as interp1d

        # Try linear interpolation; should oversample first...
        f = interp1d(t, y, bounds_error=False, fill_value=0)
        y = y(td)

    return y


def initial_value(expr, var=None):

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    return sym.limit(expr * var, var, sym.oo)


def final_value(expr, var=None):

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    return sym.limit(expr * var, var, 0)


class Zs(sExpr):
    """s-domain impedance value"""

    @classmethod
    def C(cls, Cval):
    
        return cls(1 / Cval).integrate()


    @classmethod
    def G(cls, Gval):
    
        return cls(1 / Gval)


    @classmethod
    def L(cls, Lval):
    
        return cls(Lval).differentiate()


    @classmethod
    def R(cls, Rval):
    
        return cls(Rval)


    def cpt(self):

        if self.is_number:
            return R(self.expr)

        z = self * s

        if z.is_number:
            return C((1 / z).expr)

        z = self / s

        if z.is_number:
            return L(z.expr)

        # Need a combination of components.
        return self


class Ys(sExpr):
    """s-domain admittance value"""
    

    @classmethod
    def C(cls, Cval):
    
        return cls(Cval).differentiate()


    @classmethod
    def G(cls, Gval):
    
        return cls(Gval)


    @classmethod
    def L(cls, Lval):
    
        return cls(1 / Lval).integrate()


    @classmethod
    def R(cls, Rval):
    
        return cls(1 / Rval)


    def cpt(self):

        if self.is_number:
            return G(self.expr)

        y = self * s

        if y.is_number:
            return L((1 / y).expr)

        y = self / s

        if y.is_number:
            return C(y.expr)

        # Need a combination of components.
        return self


class Vs(sExpr):
    """s-domain voltage (units V s / radian)"""
    

    def cpt(self):

        v = self * s

        if v.is_number:
            return Vdc(v.expr)

        # Need a combination of components.
        return self


class Is(sExpr):
    """s-domain current (units A s / radian)"""
    
    def cpt(self):

        i = self * s

        if i.is_number:
            return Idc(i.expr)

        # Need a combination of components.
        return self


class Avs(sExpr):
    """s-domain voltage ratio"""
    pass


class Ais(sExpr):
    """s-domain current ratio"""
    pass



class NetObject(object):

    def _tweak_args(self):

        if not hasattr(self, 'args'):
            return ()

        args = self.args
            # Drop the initial condition for L or C if it is zero.
        if isinstance(self, (L, C)) and args[1] == 0:
            args = args[:-1]
            
        modargs = []
        for arg in args:
            if isinstance(arg, sExpr):
                arg = arg.expr
            modargs.append(arg)
        return modargs


    def __repr__(self):

        argsrepr = ', '.join([arg.__repr__() for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)


    def __str__(self):
        
        return self.__repr__()


    def pretty(self):

        argsrepr = ', '.join([sym.pretty(arg) for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)


    def latex(self):

        argsrepr = ', '.join([sym.latex(arg) for arg in self._tweak_args()])
        return '\\mathrm{%s}(%s)' % (self.__class__.__name__, argsrepr)

    
    def simplify(self):

        return self


from lcapy.oneport import L, C, R, G, I, V, Idc, Vdc



