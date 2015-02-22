"""
This module provides the core functions and classes for Lcapy.

To print the rational functions in canonical form (with the highest
power of s in the denominator with a unity coefficient), use
print(x.canonical()).

For additional documentation, see the Lcapy tutorial.

Copyright 2014, 2015 Michael Hayes, UCECE
"""

from __future__ import division
import numpy as np
import sympy as sym
from sympy.utilities.lambdify import lambdify
import sys

# Note imports at bottom to avoid circular dependencies


__all__ = ('pprint', 'pretty', 'latex', 'DeltaWye', 'WyeDelta', 'tf',
           'zp2tf', 'Expr', 's', 'sExpr', 't', 'tExpr', 'f', 'fExpr', 'cExpr',
           'omega', 'omegaExpr', 'pi', 'cos', 'sin', 'exp', 'sqrt',
           'H', 'Heaviside', 'DiracDelta', 'j', 'u', 'delta',
           'Vector', 'Matrix', 'VsVector', 'IsVector', 'YsVector', 'ZsVector')


def sympify(arg, real=False):

    if isinstance(arg, (sym.symbol.Symbol, sym.symbol.Expr)):
        return arg

    if (real and isinstance(arg, str) and arg.isalnum() and
        not arg[0].isdigit()):
        arg = sym.symbols(arg, real=True)

    if isinstance(arg, str):
        # Perhaps have dictionary of functions and their replacements?
        arg = arg.replace('u(t', 'Heaviside(t')
        arg = arg.replace('delta(t', 'DiracDelta(t')
        
    # Why doesn't sympy do this?
    if isinstance(arg, complex):
        re = sym.sympify(arg.real, rational=True)
        im = sym.sympify(arg.imag, rational=True)
        if im == 1.0:
            arg = re + sym.I
        else:
            arg = re + sym.I * im
    else:
        arg = sym.sympify(arg, rational=True)

    return arg


ssym = sympify('s')
tsym = sympify('t', real=True)
fsym = sympify('f', real=True)
omegasym = sympify('omega', real=True)



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

    # Perhaps have lookup table for operands to determine
    # the resultant type?  For example, Vs / Vs -> Hs
    # Vs / Is -> Zs,  Is * Zs -> Vs

    def __init__(self, arg, real=False):

        if isinstance(arg, Expr):
            arg = arg.expr

        self.expr = sympify(arg, real=real)

    @property
    def val(self):
        """Return floating point value of expression if it can be evaluated,
        otherwise the expression"""

        return self.evalf()

    def __getattr__(self, attr):

        # This gets called if there is no explicit attribute attr for
        # this instance.  We call the method of the wrapped sympy
        # class and rewrap the returned value if it is a sympy Expr
        # object.

        expr = self.expr
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

        return self.expr.__str__()

    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.expr)

    def _repr_pretty_(self, p, cycle):

        p.pretty(self.expr)

    def _repr_latex_(self):

        return '$%s$' % self.latex()


    def __abs__(self):
        """Absolute value"""

        return self.__class__(abs(self.expr))

    def __neg__(self):
        """Negation"""

        return self.__class__(-self.expr)

    def __rdiv__(self, x):
        """Reverse divide"""

        x = self.__class__(x)
        return self.__class__(x.expr / self.expr)

    def __rtruediv__(self, x):
        """Reverse true divide"""

        x = self.__class__(x)
        return self.__class__(x.expr / self.expr)

    def __mul__(self, x):
        """Multiply"""

        x = self.__class__(x)
        return self.__class__(self.expr * x.expr)

    def __rmul__(self, x):
        """Reverse multiply"""

        x = self.__class__(x)
        return self.__class__(self.expr * x.expr)

    def __div__(self, x):
        """Divide"""

        x = self.__class__(x)
        return self.__class__(self.expr / x.expr)

    def __truediv__(self, x):
        """True divide"""

        x = self.__class__(x)
        return self.__class__(self.expr / x.expr)

    def __add__(self, x):
        """Add"""

        x = self.__class__(x)
        return self.__class__(self.expr + x.expr)

    def __radd__(self, x):
        """Reverse add"""

        x = self.__class__(x)
        return self.__class__(self.expr + x.expr)

    def __rsub__(self, x):
        """Reverse subtract"""

        x = self.__class__(x)
        return self.__class__(x.expr - self.expr)

    def __sub__(self, x):
        """Subtract"""

        x = self.__class__(x)
        return self.__class__(self.expr - x.expr)

    def __pow__(self, x):
        """Pow"""

        x = self.__class__(x)
        return self.__class__(self.expr ** x.expr)

    def __or__(self, x):
        """Parallel combination"""

        return self.parallel(x)

    def __eq__(self, x):
        """Equality"""

        if x is None:
            return False

        x = self.__class__(x)

        # This fails if one of the operands has the is_real attribute
        # end the other doesn't...
        return self.expr == x.expr

    def __ne__(self, x):
        """Inequality"""

        if x is None:
            return True

        x = self.__class__(x)
        return self.expr != x.expr

    def parallel(self, x):
        """Parallel combination"""

        x = self.__class__(x)
        return self.__class__(self.expr * x.expr / (self.expr + x.expr))

    def _pretty(self, *args, **kwargs):
        """Make pretty string"""

        # This works in conjunction with Printer._print
        # It is a hack to allow printing of _Matrix types
        # and its elements.
        expr = self.expr
        printer = args[0]

        return printer._print(expr)

    def _latex(self, *args, **kwargs):
        """Make latex string"""

        # This works in conjunction with LatexPrinter._print
        # It is a hack to allow printing of _Matrix types
        # and its elements.
        expr = self.expr
        printer = args[0]

        string = printer._print(expr)
        # sympy uses theta for Heaviside , use u(t) although I prefer H(t)
        string = string.replace(r'\theta\left', r'u\left')

        # Should replace _{xxx} with _{\mathrm{xxx}} if len(xxx) > 1

        return string

    def pretty(self):
        """Make pretty string"""
        return sym.pretty(self.expr)

    def prettyans(self, name):
        """Make pretty string with LHS name"""

        return sym.pretty(sym.Eq(sympify(name), self.expr))

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

        string = sym.latex(self.expr)
        # sympy uses theta for Heaviside
        return string.replace(r'\theta\left', r'u\left')

    def latexans(self, name):
        """Print latex string with LHS name"""

        expr = sym.Eq(sympify(name), self.expr)

        return sym.latex(expr)

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

    @property
    def conjugate(self):
        """Return complex conjugate"""

        return self.__class__(sym.conjugate(self.expr))


    @property
    def real(self):
        """Return real part"""

        return self.__class__(sym.re(self.expr).simplify())


    @property
    def imag(self):
        """Return imaginary part"""

        return self.__class__(sym.im(self.expr).simplify())

    
    def rationalize_denominator(self):
        """Rationalize denominator by multiplying numerator and denominator by
        complex conjugate of denominator"""

        N = self.N
        D = self.D
        Dconj = D.conjugate
        Nnew = (N * Dconj).simplify()
        #Dnew = (D * Dconj).simplify()
        Dnew = (D.real**2 + D.imag**2).simplify()

        Nnew = Nnew.real + j * Nnew.imag

        return Nnew / Dnew

    @property
    def magnitude(self):
        """Return magnitude"""

        x = self.rationalize_denominator()
        N = x.N
        Dnew = x.D
        Nnew = sqrt((N.real**2 + N.imag**2).simplify())

        return Nnew / Dnew

    @property
    def abs(self):
        """Return magnitude"""

        return self.magnitude

    @property
    def phase(self):
        """Return phase"""

        x = self.rationalize_denominator()
        N = x.N
        
        return self.__class__(sym.atan2(N.imag, N.real))

    @property
    def angle(self):
        """Return phase"""

        return self.phase

    @property
    def is_number(self):

        return self.expr.is_number

    @property
    def is_constant(self):

        return self.expr.is_constant()

    def evaluate(self, arg):
        """Evaluate expression at arg.
        arg may be a scalar, or a vector.

        Note, expressions such as exp(-alpha*t) * Heaviside(t) will
        not evaluate correctly since the exp will overflow for -t and
        produce an Inf.  When this is multiplied by 0 from the
        Heaviside function we get Nan. """

        # Perhaps should check if expr.args[1] == Heaviside('t') and not
        # evaluate if t < 0?

        func = lambdify(self.var, self.expr, ("numpy", "sympy", "math"))

        if np.isscalar(arg):
            v1 = arg
        else:
            v1 = arg[0]

        try:
            response = func(v1)

        except NameError:
            raise RuntimeError('Cannot evaluate expression')

        except AttributeError:
            raise RuntimeError(
                'Cannot evaluate expression,'
                ' probably have undefined symbols, such as Dirac delta')

        if np.isscalar(arg):
            return response

        try:
            response = np.array([complex(func(v1)) for v1 in arg])
        except TypeError:
            raise TypeError(
                'Cannot evaluate expression, probably have undefined symbols')

        return response

    def subs(self, arg):
        """Substitute arg for variable."""

        if isinstance(arg, Expr):
            cls = arg.__class__
            arg = arg.expr
        else:
            cls = self.__class__

        return cls(self.expr.subs(self.var, arg))

    def __call__(self, arg):

        return self.subs(arg)

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

    var = ssym

    def __init__(self, val):

        super(sExpr, self).__init__(val)
        self._laplace_conjugate_class = tExpr

        if self.expr.find(tsym) != set():
            raise ValueError(
                's-domain expression %s cannot depend on t' % self.expr)

    def _as_ratfun_delay(self):

        expr, var = self.expr, self.var

        F = sym.factor(expr).as_ordered_factors()

        delay = sympify(0)
        ratfun = sympify(1)
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
            raise ValueError('Expression not a product of rational function'
                             ' and exponential')

        numer, denom = ratfun.as_numer_denom()
        N = sym.Poly(numer, var)
        D = sym.Poly(denom, var)

        return N, D, delay

    def differentiate(self):
        """Differentiate (multiply by s)"""

        return self.__class__(self.expr * self.var)

    def integrate(self):
        """Integrate (divide by s)"""

        return self.__class__(self.expr / self.var)

    def delay(self, T):
        """Apply delay of T seconds by multiplying by exp(-s T)"""

        T = self.__class__(T)
        return self.__class__(self.expr * sym.exp(-s * T))

    def jomega(self):
        """Return expression with s = j omega"""

        w = omegaExpr(omegasym)
        return self.subs(sym.I * w)

    def roots(self):
        """Return roots of expression as a dictionary
        Note this may not find them all."""

        return Exprdict(sym.roots(sym.Poly(self.expr, self.var)))

    def zeros(self):
        """Return zeroes of expression as a dictionary
        Note this may not find them all."""

        return self.N.roots()

    def poles(self):
        """Return poles of expression as a dictionary
        Note this may not find them all."""

        return self.D.roots()

    def residue(self, pole, poles):

        expr, var = self.expr, self.var

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

    def _as_residue_parts(self):
        """Return residues of expression"""

        var = self.var
        N, D, delay = self._as_ratfun_delay()

        Q, M = N.div(D)

        expr = M / D
        sexpr = sExpr(expr)

        P = sexpr.poles()
        F = []
        R = []
        for p in P:

            # Number of occurrences of the pole.
            N = P[p]

            f = var - p

            if N == 1:
                F.append(f)
                R.append(sexpr.residue(p, P))
                continue

            # Handle repeated poles.
            expr2 = expr * f ** N
            for n in range(1, N + 1):
                m = N - n
                F.append(f ** n)
                dexpr = sym.diff(expr2, var, m)
                R.append(sym.limit(dexpr, var, p) / sym.factorial(m))

        return F, R, Q, delay

    def canonical(self):
        """Convert rational function to canonical form with unity
        highest power of denominator.

        See also general, partfrac, mixedfrac, and ZPK"""

        N, D, delay = self._as_ratfun_delay()

        K = sym.cancel(N.LC() / D.LC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        # Divide by leading coefficient
        N = N.monic()
        D = D.monic()

        expr = K * (N / D)

        return self.__class__(expr)

    def general(self):
        """Convert rational function to general form.

        See also canonical, partfrac, mixedfrac, and ZPK"""

        N, D, delay = self._as_ratfun_delay()

        expr = sym.cancel(N / D, self.var)
        if delay != 0:
            expr *= sym.exp(self.var * delay)

        return self.__class__(expr)

    def partfrac(self):
        """Convert rational function into partial fraction form.

        See also canonical, mixedfrac, general, and ZPK"""

        F, R, Q, delay = self._as_residue_parts()

        expr = Q.as_expr()
        for f, r in zip(F, R):
            expr += r / f

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return self.__class__(expr)

    def mixedfrac(self):
        """Convert rational function into mixed fraction form.

        See also canonical, general, partfrac and ZPK"""

        N, D, delay = self._as_ratfun_delay()

        Q, M = N.div(D)

        expr = Q + M / D

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return self.__class__(expr)

    def ZPK(self):
        """Convert to pole-zero-gain (PZK) form.

        See also canonical, general, mixedfrac, and partfrac"""

        N, D, delay = self._as_ratfun_delay()

        K = sym.cancel(N.LC() / D.LC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        zeros = sym.roots(N)
        poles = sym.roots(D)

        return self.__class__(_zp2tf(zeros, poles, K, self.var))

    def initial_value(self):
        """Determine value at t = 0"""

        return self.__class__(sym.limit(self.expr * self.var, self.var, sym.oo))

    def final_value(self):
        """Determine value at t = oo"""

        return self.__class__(sym.limit(self.expr * self.var, self.var, 0))


    def _inverse_laplace(self):

        var = self.var
        N, D, delay = self._as_ratfun_delay()

        Q, M = N.div(D)

        result1 = 0

        # Delayed time.
        td = tsym - delay

        if Q:
            C = Q.all_coeffs()
            for n, c in enumerate(C):
                result1 += c * sym.diff(sym.DiracDelta(td), tsym, len(C) - n - 1)

        expr = M / D
        sexpr = sExpr(expr)

        P = sexpr.poles()
        result2 = 0

        P2 = P.copy()

        for p in P2:

            # Number of occurrences of the pole.
            N = P2[p]

            if N == 0:
                continue

            f = var - p

            if N == 1:
                r = sexpr.residue(p, P)

                pc = p.conjugate()
                if pc != p and pc in P:
                    # Remove conjugate from poles and process pole with its
                    # conjugate
                    P2[pc] = 0
                  
                    p_re = sym.re(p)
                    p_im = sym.im(p)
                    r_re = sym.re(r)
                    r_im = sym.im(r)
                    etd = sym.exp(p_re * td)
                    result2 += 2 * r_re * etd * sym.cos(p_im * td)
                    result2 += 2 * r_im * etd * sym.sin(p_im * td)
                else:
                    result2 += r * sym.exp(p * td)
                continue

            # Handle repeated poles.
            expr2 = expr * f ** N
            for n in range(1, N + 1):
                m = N - n
                r = sym.limit(
                    sym.diff(expr2, var, m), var, p) / sym.factorial(m)
                result2 += r * sym.exp(p * td) * td**(n - 1)

        return result1 + result2 * sym.Heaviside(td)

    def inverse_laplace(self):
        """Attempt inverse Laplace transform"""

        try:
            result = self._inverse_laplace()
        except:

            print('Determining inverse Laplace transform with sympy...')

            # Try splitting into partial fractions to help sympy.
            expr = self.partfrac().expr

            # This barfs when needing to generate Dirac deltas
            from sympy.integrals.transforms import inverse_laplace_transform
            result = inverse_laplace_transform(expr, t, self.var)

        if hasattr(self, '_laplace_conjugate_class'):
            result = self._laplace_conjugate_class(result)
        return result

    def transient_response(self, t=None):
        """Evaluate transient (impulse) response"""

        texpr = self.inverse_laplace()

        if t is None:
            return texpr

        print('Evaluating inverse Laplace transform...')
        return texpr.evaluate(t)

    def impulse_response(self, t=None):
        """Evaluate transient (impulse) response"""

        return self.transient_response(t)

    def step_response(self, t=None):
        """Evaluate step response"""

        return (self / self.var).transient_response(t)

    def frequency_response(self, f=None):
        """Convert to frequency domain and evaluate response is frequency
        specified"""

        if f is None:
            return self.subs(sym.I * 2 * sym.pi * fsym)

        return self.evaluate(2j * np.pi * f)

    def response(self, x, t):
        """Evaluate response to input signal x at times t"""

        N, D, delay = self._as_ratfun_delay()

        Q, M = N.div(D)
        expr = M / D

        h = sExpr(expr).transient_response(t)
        y = np.convolve(x, h)

        if Q:
            # Handle Dirac Deltas and derivatives.
            C = Q.all_coeffs()
            dt = np.diff(t)
            for n, c in enumerate(C):

                y += c * x

                x = np.diff(x) / dt
                x = np.hstack((x, 0))

        from scipy.interpolate import interp1d

        # Try linear interpolation; should oversample first...
        y = interp1d(t, y, bounds_error=False, fill_value=0)
        td = t - delay
        y = y(td)

        return y

    def decompose(self):

        N, D, delay = self._as_ratfun_delay()

        return N, D, delay

    def evaluate(self, svector):

        return super(sExpr, self).evaluate(svector)

    def plot(self, t=None, **kwargs):

        from lcapy.plot import plot_pole_zero

        plot_pole_zero(self, **kwargs)


class fExpr(Expr):

    """Fourier domain expression or symbol"""

    var = fsym
    domain_name = 'Frequency'
    domain_units = 'Hz'

    def __init__(self, val, real=False):

        super(fExpr, self).__init__(val, real=real)
        self._fourier_conjugate_class = tExpr

        if self.expr.find(ssym) != set():
            raise ValueError(
                'f-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find(tsym) != set():
            raise ValueError(
                'f-domain expression %s cannot depend on t' % self.expr)

    def inverse_fourier(self):
        """Attempt inverse Fourier transform"""

        result = sym.inverse_fourier_transform(self.expr, tsym, self.var.expr)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        return result

    def plot(self, f=None, **kwargs):

        from lcapy.plot import plot_frequency
        plot_frequency(self, f, **kwargs)


class omegaExpr(Expr):

    """Fourier domain expression or symbol (angular frequency)"""

    var = omegasym
    domain_name = 'Angular frequency'
    domain_units = 'rad/s'

    def __init__(self, val, real=False):

        super(omegaExpr, self).__init__(val, real=real)
        self._fourier_conjugate_class = tExpr

        if self.expr.find(ssym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find(tsym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on t' % self.expr)

    def inverse_fourier(self):
        """Attempt inverse Fourier transform"""

        result = sym.inverse_fourier_transform(self.expr, tsym, self.var.expr)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        return result

    def plot(self, f=None, **kwargs):

        from lcapy.plot import plot_angular_frequency
        plot_angular_frequency(self, f, **kwargs)


class tExpr(Expr):

    """t-domain expression or symbol"""

    var = tsym
    domain_name = 'Time'
    domain_units = 's'

    def __init__(self, val, real=True):

        super(tExpr, self).__init__(val, real=real)
        self._fourier_conjugate_class = fExpr
        self._laplace_conjugate_class = sExpr

        if self.expr.find(ssym) != set():
            raise ValueError(
                't-domain expression %s cannot depend on s' % self.expr)

    def laplace(self):
        """Attempt Laplace transform"""
        
        var = self.var
        if var == tsym:
            # Hack since sympy gives up on Laplace transform if t real!
            var = sym.symbols('t')

        F, a, cond = sym.laplace_transform(self, var, ssym)

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

        response = super(tExpr, self).evaluate(tvector)
        if np.iscomplexobj(response) and np.allclose(response.imag, 0.0):
            response = response.real
        return response

    def plot(self, t=None, **kwargs):

        from lcapy.plot import plot_time
        plot_time(self, t, **kwargs)


class cExpr(Expr):

    """Constant real expression or symbol"""

    def __init__(self, val):

        # FIXME for expressions of cExpr
        if False and not isinstance(val, (cExpr, int, float, str)):
            raise ValueError(
                '%s of type %s not int, float, or str' % (val, type(val)))

        super(cExpr, self).__init__(val, real=True)


s = sExpr('s')
t = tExpr('t', real=True)
f = fExpr('f', real=True)
omega = omegaExpr('omega', real=True)
pi = sym.pi
j = sym.I


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
    # original type information and associated methods.  As a hack, we
    # try to wrap elements when they are read using __getitem__.  This
    # assumes that all the elements have the same type.  This is not
    # the case for A, B, G, and H matrices.  This could he handled by
    # having another matrix to specify the type for each element.

    _typewrap = sExpr

    def __getitem__(self, key):

        item = super(Matrix, self).__getitem__(key)

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

        from copy import copy
        new = copy(self)

        for i in range(self.rows):
            for j in range(self.cols):
                new[i, j] = getattr(self[i, j], method)()

        return new

    def canonical(self):

        return self._reformat('canonical')

    def general(self):

        return self._reformat('general')

    def mixedfrac(self):

        return self._reformat('mixedfrac')

    def partfrac(self):

        return self._reformat('partfrac')

    def ZPK(self):

        return self._reformat('ZPK')


class Vector(Matrix):

    def __new__(cls, *args):

        args = [sympify(arg) for arg in args]

        if len(args) == 2:
            return super(Vector, cls).__new__(cls, (args[0], args[1]))

        return super(Vector, cls).__new__(cls, *args)


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
        var = ssym

    N = sym.Poly(numer, var)
    D = sym.Poly(denom, var)

    return Hs(N / D)


def _zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists of zeros and poles,
    and from a constant gain"""

    if var is None:
        var = ssym

    zeros = sympify(zeros)
    poles = sympify(poles)

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

    return Hs(_zp2tf(zeros, poles, K, var))


# Perhaps use a factory to create the following classes?

class Zs(sExpr):

    """s-domain impedance value"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val):

        super(Zs, self).__init__(val)
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

        super(Ys, self).__init__(val)
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

        super(Vs, self).__init__(val)
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

        super(Is, self).__init__(val)
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

        super(Hs, self).__init__(val)
        self._laplace_conjugate_class = Ht


class Yt(tExpr):

    """t-domain 'admittance' value"""

    units = 'siemens/s'

    def __init__(self, val):

        super(Yt, self).__init__(val)
        self._laplace_conjugate_class = Ys
        self._fourier_conjugate_class = Yf


class Zt(tExpr):

    """t-domain 'impedance' value"""

    units = 'ohms/s'

    def __init__(self, val):

        super(Zt, self).__init__(val)
        self._laplace_conjugate_class = Zs
        self._fourier_conjugate_class = Zf


class Vt(tExpr):

    """t-domain voltage (units V)"""

    quantity = 'Voltage'
    units = 'V'

    def __init__(self, val):

        super(Vt, self).__init__(val)
        self._laplace_conjugate_class = Vs
        self._fourier_conjugate_class = Vf


class It(tExpr):

    """t-domain current (units A)"""

    quantity = 'Current'
    units = 'A'

    def __init__(self, val):

        super(It, self).__init__(val)
        self._laplace_conjugate_class = Is
        self._fourier_conjugate_class = If


class Ht(tExpr):

    """impulse response"""

    quantity = 'Impulse response'
    units = '1/s'

    def __init__(self, val):

        super(Ht, self).__init__(val)
        self._laplace_conjugate_class = Hs
        self._fourier_conjugate_class = Hf


class Yf(fExpr):

    """f-domain admittance"""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val):

        super(Yf, self).__init__(val)
        self._fourier_conjugate_class = Yt


class Zf(fExpr):

    """f-domain impedance"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val):

        super(Zf, self).__init__(val)
        self._fourier_conjugate_class = Zt


class Vf(fExpr):

    """f-domain voltage (units V)"""

    quantity = 'Voltage spectrum'
    units = 'V/Hz'

    def __init__(self, val):

        super(Vf, self).__init__(val)
        self._fourier_conjugate_class = Vt


class If(fExpr):

    """f-domain current (units A)"""

    quantity = 'Current spectrum'
    units = 'A/Hz'

    def __init__(self, val):

        super(If, self).__init__(val)
        self._fourier_conjugate_class = It


class Hf(fExpr):

    """d-domain transfer function response"""

    quantity = 'Transfer function'
    units = '1'

    def __init__(self, val):

        super(Hf, self).__init__(val)
        self._fourier_conjugate_class = Ht


class VsVector(Vector):

    _typewrap = Vs


class IsVector(Vector):

    _typewrap = Is


class YsVector(Vector):

    _typewrap = Ys


class ZsVector(Vector):

    _typewrap = Zs


class NetObject(object):

    def __init__(self, args):

        self.args = args

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


def _funcwrap(func, *args):

    expr = args[0]
    cls = expr.__class__

    if isinstance(expr, Expr):
        expr = expr.expr

    return cls(func(expr, *args[1:]))


def sin(expr):

    return _funcwrap(sym.sin, expr)


def cos(expr):

    return _funcwrap(sym.cos, expr)


def exp(expr):

    return _funcwrap(sym.exp, expr)


def sqrt(expr):

    return _funcwrap(sym.sqrt, expr)


def Heaviside(expr):
    """Heaviside's unit step"""

    return _funcwrap(sym.Heaviside, expr)


def H(expr):
    """Heaviside's unit step"""

    return Heaviside(expr)


def u(expr):
    """Heaviside's unit step"""

    return Heaviside(expr)


def DiracDelta(*args):
    """Dirac delta (impulse)"""

    return _funcwrap(sym.DiracDelta, *args)


def delta(expr, *args):
    """Dirac delta (impulse)"""

    return DiracDelta(expr, *args)


from lcapy.oneport import L, C, R, G, Idc, Vdc
