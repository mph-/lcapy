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
import sys

# Note imports at bottom to avoid circuit dependencies


__all__ = ('pprint', 'pretty', 'latex', 'DeltaWye', 'WyeDelta', 'tf', 
           'zp2tf', 'poles', 'zeros', 'residue', 'residues', 'partfrac',
           'general', 'canonical', 'ZPK', 'inverse_laplace', 'initial_value',
           'transient_response', 'response', 'final_value', 'Expr', 
           's', 'sExpr',  't', 'tExpr', 'f', 'fExpr', 'cExpr',
           'omega', 'omegaExpr', 'pi', 'cos', 'sin', 'exp', 'sqrt',
           'H', 'DiracDelta',
           'Vector', 'Matrix', 'VsVector', 'IsVector', 'YVector', 'ZVector')


class Exprdict(dict):
    """Decorator class for dictionary created by sympy"""

    def pprint(self):
        """Pretty print"""
        
        return sym.pprint(self)


    def latex(self):
        """Latex"""

        return sym.latex(self)


class Expr(object):
    """Decorator class for sympy classes derived from sympy.Expr"""
    

    @property
    def expr(self):    
        return self.val
    
    
    def __init__(self, val, real=False):
        
        if isinstance(val, sExpr):
            val = val.val
            
        if real and isinstance(val, str) and val.isalnum() and not val[0].isdigit():
            val = sym.symbols(val, real=True)

        if isinstance(val, str):
            val = val.replace('u(t', 'Heaviside(t')
            val = val.replace('delta(t', 'DiracDelta(t')
            
        val = sym.sympify(val, rational=True)

        self.val = val


    def __getattr__(self, attr):

        # This gets called if there is no explicit attribute attr for
        # this instance.  We call the method of the wrapped sympy
        # class and rewrap the returned value if it is a sympy Expr
        # object.

        expr = self.val
        if hasattr(expr, attr):
            
            def wrap(*args):

                ret = getattr(expr, attr)(*args)
                if not isinstance(ret, sym.Expr):
                    return ret

                # Wrap the return value
                return self.__class__(ret)

            return wrap

            
        # Try looking for a sympy function with the same name,
        # such as sqrt, log, etc.
        # On second thoughts, this may confuse the user since
        # we will pick up methods such as laplace_transform.
        # Perhaps should have a list of allowable functions?
        if True or not hasattr(sym, attr):
            raise AttributeError(
                "%s has no attribute %s." % (self.__class__.__name__, attr))
            
        def wrap1(*args):

            ret = getattr(sym, attr)(expr, *args)
            if not isinstance(ret, sym.Expr):
                return ret

            # Wrap the return value
            return self.__class__(ret)

        return wrap1


    def __str__(self):

        return self.val.__str__()


    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.val)


    def _repr_pretty_(self, p, cycle):

        p.pretty(self.expr)


    def _repr_latex_(self):

        return '$%s$' % self.latex()

    
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

        if x is None:
            return False

        x = self.__class__(x)
        return self.val == x.val


    def __ne__(self, x):
        """Inequality"""

        if x is None:
            return True

        x = self.__class__(x)
        return self.val != x.val


    def parallel(self, x):
        """Parallel combination"""
        
        x = self.__class__(x)
        return self.__class__(self.val * x.val / (self.val + x.val))
    

    def _pretty(self, *args, **kwargs):
        """Make pretty string"""
        
        # This works in conjunction with Printer._print
        # It is a hack to allow printing of _Matrix types 
        # and its elements.
        expr = self.val
        printer = args[0]

        return printer._print(expr)


    def _latex(self, *args, **kwargs):
        """Make latex string"""
        
        # This works in conjunction with LatexPrinter._print
        # It is a hack to allow printing of _Matrix types 
        # and its elements.
        expr = self.val
        printer = args[0]

        string = printer._print(expr)
        # sympy uses theta for Heaviside , use u(t) although I prefer H(t)
        string = string.replace(r'\theta\left', r'u\left')

        # Should replace _{xxx} with _{\mathrm{xxx}} if len(xxx) > 1

        return string



    def pretty(self):
        """Make pretty string"""
        return sym.pretty(self.val)


    def prettyans(self, name):
        """Make pretty string with LHS name"""

        return sym.pretty(sym.Eq(sym.sympify(name), self.val))


    def pprint(self):
        """Pretty print"""
        
        # If interactive use pretty, otherwise use latex
        if hasattr(sys, 'ps1'):
            print(self.pretty())
        else:
            print(self.latex())


    def pprintans(self, name):
        """Pretty print string with LHS name"""
        print(self.prettyans(name))


    def latex(self):
        """Make latex string"""

        string = sym.latex(self.val)
        # sympy uses theta for Heaviside 
        return string.replace(r'\theta\left', r'u\left')


    def latexans(self, name):
        """Print latex string with LHS name"""

        expr = sym.Eq(sym.sympify(name), self.val)

        return sym.latex(expr)


    @property
    def is_number(self):

        return self.expr.is_number


    @property
    def is_constant(self):

        return self.expr.is_constant()


    def evaluate(self, vector, var):

        func = lambdify(var, self.expr, ("numpy", "sympy", "math"))
        
        # Expressions such as exp(-alpha*t) * Heaviside(t)
        # will not evaluate correctly since the exp will overflow
        # for -t and produce an Inf.  When this is multiplied by
        # 0 from the Heaviside function we get Nan.
        # Perhaps should check if expr.args[1] == Heaviside('t')
        # and not evaluate if t < 0?

        if np.isscalar(vector):
            v1 = vector
        else:
            v1 = vector[0]


        try:
            response = func(v1)

        except NameError:
            raise RuntimeError('Cannot evaluate expression')

        except AttributeError:
            raise RuntimeError('Cannot evaluate expression, probably have undefined symbols, such as Dirac delta')

        if np.isscalar(vector):
            return response

        try:
            response = np.array([complex(func(v1)) for v1 in vector])
        except TypeError:
            raise TypeError('Cannot evaluate expression, probably have undefined symbols')
            
        return response


    def __call__(self, vector):

        return self.evaluate(vector)


    @property
    def label(self):

        label = ''
        if hasattr(self, 'quantity'):
            label += '%s' % self.quantity
        if hasattr(self, 'units'):
            label += ' (%s)' % self.units
        return label


    @property
    def domain_label(self):

        label = ''
        if hasattr(self, 'domain_name'):
            label += '%s' % self.domain_name
        if hasattr(self, 'domain_units'):
            label += ' (%s)' % self.domain_units
        return label


class sExpr(Expr):
    """s-domain expression or symbol"""
    
    var = sym.symbols('s')

    def __init__(self, val):

        super (sExpr, self).__init__(val)
        self._laplace_conjugate_class = tExpr
        
        if self.expr.find('t') != set():
            raise ValueError('s-domain expression %s cannot depend on t' % self.expr)


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


    def omega(self):
        """Return expression with s = j omega"""

        omega = sym.symbols('omega')
        return omegaExpr(self.subs(s, sym.I * omega))


    def zeros(self):
        """Return zeroes of expression as a dictionary"""
        
        return Exprdict(zeros(self.expr, self.var))


    def poles(self):
        """Return poles of expression as a dictionary"""
        
        return Exprdict(poles(self.expr, self.var))


    def residues(self):
        """Return residues of expression"""
        
        return residues(self.expr, self.var)


    def canonical(self):
        """Convert rational function to canonical form with unity
        highest power of denominator.

        See also general, partfrac, mixedfrac, and ZPK"""
        
        return self.__class__(canonical(self.expr, self.var))


    def general(self):
        """Convert rational function to general form.

        See also canonical, partfrac, mixedfrac, and ZPK"""
        
        return self.__class__(general(self.expr, self.var))


    def partfrac(self):
        """Convert rational function into partial fraction form.

        See also canonical, mixedfrac, general, and ZPK"""

        return self.__class__(partfrac(self.expr, self.var))


    def mixedfrac(self):
        """Convert rational function into mixed fraction form.

        See also canonical, general, partfrac and ZPK"""

        return self.__class__(mixedfrac(self.expr, self.var))


    def ZPK(self):
        """Convert to pole-zero-gain (PZK) form.
        
        See also canonical, general, mixedfrac, and partfrac"""
        
        return self.__class__(ZPK(self.expr, self.var))


    def initial_value(self):
        """Determine value at t = 0"""
        
        return initial_value(self.expr, self.var)


    def final_value(self):
        """Determine value at t = oo"""
        
        return final_value(self.expr, self.var)


    @property
    def N(self):
        """Return numerator of rational function"""

        expr = self.expr
        if not expr.is_rational_function(self):
            raise ValueError('Expression not a rational function')

        numer, denom = expr.as_numer_denom()
        return self.__class__(numer)


    @property
    def D(self):
        """Return denominator of rational function"""

        expr = self.expr
        if not expr.is_rational_function(self):
            raise ValueError('Expression not a rational function')

        numer, denom = expr.as_numer_denom()
        return self.__class__(denom)


    def _as_ratfun_parts(self):
        
        return _as_ratfun_parts(self.expr, self.var)


    def inverse_laplace(self):
        """Attempt inverse Laplace transform"""
        
        result = inverse_laplace(self.expr, t, self.var)
        if hasattr(self, '_laplace_conjugate_class'):
            result = self._laplace_conjugate_class(result)
        return result


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
        
        if f is None:
            fsym = sym.symbols('f')
            return fExpr(self.subs(s, sym.I * 2 * sym.pi * fsym))

        return self.evaluate(2j * np.pi * f)


    def response(self, x, t):
        """Evaluate response to input signal x"""

        return response(self.expr, x, t, self.var)


    def decompose(self):

        ratfun, delay = _as_ratfun_delay(self.expr, self.var)
        
        N, D = _as_ratfun_parts(ratfun, self.var)

        return N, D, delay


    def evaluate(self, svector):

        return super (sExpr, self).evaluate(svector, sym.symbols('s'))


class fExpr(Expr):
    """Fourier domain expression or symbol"""
    
    var = sym.symbols('f')
    domain_name = 'Frequency'
    domain_units = 'Hz'

    def __init__(self, val):

        super (fExpr, self).__init__(val)
        self._fourier_conjugate_class = tExpr
        
        if self.expr.find('s') != set():
            raise ValueError('f-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find('t') != set():
            raise ValueError('f-domain expression %s cannot depend on t' % self.expr)


    def inverse_fourier(self):
        """Attempt inverse Fourier transform"""
        
        result = sym.inverse_fourier_transform(self.expr, t, self.var)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        return result


class omegaExpr(Expr):
    """Fourier domain expression or symbol (angular frequency)"""
    
    var = sym.symbols('omega')
    domain_name = 'Angular frequency'
    domain_units = 'rad/s'

    def __init__(self, val):

        super (omegaExpr, self).__init__(val)
        self._fourier_conjugate_class = tExpr
        
        if self.expr.find('s') != set():
            raise ValueError('omega-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find('t') != set():
            raise ValueError('omega-domain expression %s cannot depend on t' % self.expr)


    def inverse_fourier(self):
        """Attempt inverse Fourier transform"""
        
        result = sym.inverse_fourier_transform(self.expr, t, self.var)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        return result


class tExpr(Expr):
    """t-domain expression or symbol"""

    var = sym.symbols('t')
    domain_name = 'Time'
    domain_units = 's'

    def __init__(self, val):

        super (tExpr, self).__init__(val)
        self._fourier_conjugate_class = fExpr
        self._laplace_conjugate_class = sExpr
        
        if self.expr.find('s') != set():
            raise ValueError('t-domain expression %s cannot depend on s' % self.expr)


    def laplace(self):
        """Attempt Laplace transform"""
        
        F, a, cond = sym.laplace_transform(self.expr, self.var, s)

        if hasattr(self, '_laplace_conjugate_class'):
            F = self._laplace_conjugate_class(F)
        return F


    def fourier(self):
        """Attempt Fourier transform"""
        
        F = sym.fourier_transform(self.expr, self.var, f)

        if hasattr(self, '_fourier_conjugate_class'):
            F = self._fourier_conjugate_class(F)
        return F


    def evaluate(self, tvector):

        response = super (tExpr, self).evaluate(tvector, sym.symbols('t'))
        if np.iscomplexobj(response) and np.allclose(response.imag, 0.0):
            response = response.real
        return response


    def plot(self):

        from matplotlib.pyplot import figure
        
        # FIXME, determine useful time range...
        t = np.linspace(-0.2, 2, 400)
        v = self(t)

        fig = figure()
        ax = fig.add_subplot(111)
        ax.plot(t, v)
        ax.set_xlabel(self.domain_label)
        ax.set_ylabel(self.label)
        ax.grid(True)


class cExpr(Expr):
    """Constant expression or symbol"""

    def __init__(self, val):
        
        # FIXME for expressions of cExpr
        if False and not isinstance(val, (cExpr, int, float, str)):
            raise ValueError('%s of type %s not int, float, or str' % (val, type(val)))

        super (cExpr, self).__init__(val, real=True)


s = sExpr('s')
t = tExpr('t')
f = fExpr('f')
omega = omegaExpr('omega')
pi = sym.pi

def _guess_var(expr, var):

    if hasattr(expr, 'expr'):
        return expr.expr, expr.var

    if var is not None:
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


class Matrix(sym.Matrix):

    # Unlike numpy.ndarray, the sympy.Matrix runs all the elements
    # through sympify, creating sympy objects and thus losing the
    # original type information and associated methods.
    # As a hack, we try to wrap elements when they are read
    # using __getitem__.

    _typewrap = sExpr


    def __getitem__(self, key):

        item = super (Matrix, self).__getitem__(key)

        # The following line is to handle slicing used
        # by latex method.
        if isinstance(item, sym.Matrix):
            return item

        if hasattr(self, '_typewrap'):
            return self._typewrap(item)

        return item


    def pprint(self):

        return sym.pprint(self)


    def latex(self):

        return sym.latex(self)


    def _reformat(self, method):
        """Helper method for reformatting expression"""

        new = copy(self)

        for i in range(self.rows):
            for j in range(self.cols):
                new[i, j] = getattr(self[i, j], method)()

        return new


    def canonical(self):

        return self._reformat('canonical');


    def general(self):

        return self._reformat('general');


    def mixedfrac(self):

        return self._reformat('mixedfrac');


    def partfrac(self):

        return self._reformat('partfrac');


    def ZPK(self):

        return self._reformat('ZPK');


class Vector(Matrix):

    def __new__ (cls, *args):

        args = [sym.sympify(arg) for arg in args]

        if len(args) == 2:
            return super (Vector, cls).__new__(cls, (args[0], args[1]))

        return super (Vector, cls).__new__(cls, *args)


def DeltaWye(Z1, Z2, Z3):

    ZZ = (Z1 * Z2 + Z2 * Z3 + Z3 * Z1)
    return (ZZ / Z1, ZZ / Z2, ZZ / Z3)


def WyeDelta(Z1, Z2, Z3):

    ZZ = Z1 + Z2 + Z3
    return (Z2 * Z3 / ZZ, Z1 * Z3 / ZZ, Z1 * Z2 / ZZ)


def tf(numer, denom=1, var=None):
    """Create a transfer function from lists of the coefficient
    for the numerator and denominator"""

    if var is None:
        var = sym.symbols('s')

    N = sym.Poly(numer, var)
    D = sym.Poly(denom, var)

    return sExpr(N / D)


def _zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists of zeros and poles, 
    and from a constant gain"""
    
    if var is None:
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


def mixedfrac(expr, var=None):
    """Convert rational function into mixed fraction form.

    See also canonical, general, partfrac, and ZPK"""

    expr, var = _guess_var(expr, var)

    ratfun, delay = _as_ratfun_delay(expr, var)

    N, D = _as_ratfun_parts(ratfun, var)

    Q, M = N.div(D)

    expr = Q + M / D

    if delay != 0:
        expr *= sym.exp(-var * delay)

    return expr


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

    if not expr.is_rational_function(var):
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
            
    if not ratfun.is_rational_function(var):
        raise ValueError('Expression not a product of rational function and exponential')

    return ratfun, delay


def general(expr, var=None):
    """Convert rational function into general form.

    See also canonical, partfrac, mixedfrac, and ZPK"""

    expr, var = _guess_var(expr, var)

    ratfun, delay = _as_ratfun_delay(expr, var)

    result = sym.cancel(ratfun, var)
    if delay != 0:
        result = result * sym.exp(var * delay)

    return result


def canonical(expr, var=None):
    """Convert rational function into canonical form.

    See also general, partfrac, mixedfrac, and ZPK"""

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

    See also canonical, general, partfrac, and mixedfrac"""

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

    P2 = P.copy()

    for p in P2:
        
        # Number of occurrences of the pole.
        N = P2[p]

        if N == 0:
            continue

        f = var - p

        if N == 1:
            r = residue(expr, var, p, P)

            pc = p.conjugate()
            if pc != p and P.has_key(pc):
                # Remove conjugate from poles and process pole with its conjugate
                P2[pc] = 0
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
    if t is None:
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
    if t is None:
        return texpr

    print('Evaluating inverse Laplace transform...')
        
    return texpr.evaluate(t)


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


# Perhaps use a factory to create the following classes?

class Zs(sExpr):
    """s-domain impedance value"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val):

        super (Zs, self).__init__(val)
        self._laplace_conjugate_class = Zt


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
    
    quantity = 'Admittance'
    units = 'siemens'


    def __init__(self, val):

        super (Ys, self).__init__(val)
        self._laplace_conjugate_class = Yt


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

    quantity = 's-Voltage'
    units = 'V/Hz'    

    def __init__(self, val):

        super (Vs, self).__init__(val)
        self._laplace_conjugate_class = Vt


    def cpt(self):

        v = self * s

        if v.is_number:
            return Vdc(v.expr)

        # Need a combination of components.
        return self


class Is(sExpr):
    """s-domain current (units A s / radian)"""

    quantity = 's-Current'
    units = 'A/Hz'    


    def __init__(self, val):

        super (Is, self).__init__(val)
        self._laplace_conjugate_class = It

    
    def cpt(self):

        i = self * s

        if i.is_number:
            return Idc(i.expr)

        # Need a combination of components.
        return self


class Hs(sExpr):
    """s-domain ratio"""
    pass

    quantity = 's-ratio'
    units = ''    

    def __init__(self, val):

        super (Hs, self).__init__(val)
        self._laplace_conjugate_class = Ht


class Yt(tExpr):
    """t-domain 'admittance' value"""

    units = 'siemens/s'

    def __init__(self, val):

        super (Yt, self).__init__(val)
        self._laplace_conjugate_class = Ys
        self._fourier_conjugate_class = Yf


class Zt(tExpr):
    """t-domain 'impedance' value"""

    units = 'ohms/s'

    def __init__(self, val):

        super (Zt, self).__init__(val)
        self._laplace_conjugate_class = Zs
        self._fourier_conjugate_class = Zf


class Vt(tExpr):
    """t-domain voltage (units V)"""

    quantity = 'Voltage'
    units = 'V'

    def __init__(self, val):

        super (Vt, self).__init__(val)
        self._laplace_conjugate_class = Vs
        self._fourier_conjugate_class = Vf


class It(tExpr):
    """t-domain current (units A)"""
    
    quantity = 'Current'
    units = 'A'

    def __init__(self, val):

        super (It, self).__init__(val)
        self._laplace_conjugate_class = Is
        self._fourier_conjugate_class = If


class Ht(tExpr):
    """impulse response"""
    
    quantity = 'Impulse response'
    units = '1/s'

    def __init__(self, val):

        super (Ht, self).__init__(val)
        self._laplace_conjugate_class = Hs
        self._fourier_conjugate_class = Hf


class Yf(fExpr):
    """f-domain admittance"""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val):

        super (Yf, self).__init__(val)
        self._fourier_conjugate_class = Yt


class Zf(fExpr):
    """f-domain impedance"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val):

        super (Zf, self).__init__(val)
        self._fourier_conjugate_class = Zt


class Vf(fExpr):
    """f-domain voltage (units V)"""

    quantity = 'Voltage spectrum'
    units = 'V/Hz'

    def __init__(self, val):

        super (Vf, self).__init__(val)
        self._fourier_conjugate_class = Vt


class If(fExpr):
    """f-domain current (units A)"""
    
    quantity = 'Current spectrum'
    units = 'A/Hz'

    def __init__(self, val):

        super (If, self).__init__(val)
        self._fourier_conjugate_class = It


class Hf(fExpr):
    """d-domain transfer function response"""
    
    quantity = 'Transfer function'
    units = '1'

    def __init__(self, val):

        super (Ht, self).__init__(val)
        self._fourier_conjugate_class = Ht



class VsVector(Vector):
    
    _typewrap = Vs


class IsVector(Vector):
    
    _typewrap = Is


class YVector(Vector):

    _typewrap = Ys


class ZVector(Vector):

    _typewrap = Zs


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

        def fmt(arg):
            if False and isinstance(arg, str):
                return "'" + arg + "'" 
            return arg.__str__()
        
        argsrepr = ', '.join([fmt(arg) for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)


    def _repr_pretty_(self, p, cycle):

        p.text(self.pretty())


    def _repr_latex_(self):

        return '$%s$' % self.latex()


    def pretty(self):

        argsrepr = ', '.join([sym.pretty(arg) for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)


    def latex(self):

        argsrepr = ', '.join([sym.latex(arg) for arg in self._tweak_args()])
        return '\\mathrm{%s}(%s)' % (self.__class__.__name__, argsrepr)

    
    def simplify(self):

        return self


def sin(expr):
    
    return expr.__class__(sym.sin(expr))


def cos(expr):
    
    return expr.__class__(sym.cos(expr))


def exp(expr):
    
    return expr.__class__(sym.exp(expr))


def sqrt(expr):
    
    return expr.__class__(sym.sqrt(expr))


def H(expr):
    
    return expr.__class__(sym.Heaviside(expr))


def DiracDelta(expr):
    
    return expr.__class__(sym.DiracDelta(expr))



from lcapy.oneport import L, C, R, G, I, V, Idc, Vdc
