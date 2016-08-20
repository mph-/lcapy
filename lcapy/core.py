"""
This module provides the core functions and classes for Lcapy.

To print the rational functions in canonical form (with the highest
power of s in the denominator with a unity coefficient), use
print(x.canonical()).

For additional documentation, see the Lcapy tutorial.

Copyright 2014, 2015, 2016 Michael Hayes, UCECE
"""

from __future__ import division
from lcapy.latex import latex_str
from lcapy.acdc import is_dc, is_ac, is_causal, ACChecker
from lcapy.sympify import canonical_name, sympify1, symbols_find
import numpy as np
from sympy.core.mul import _unevaluated_Mul as uMul
from sympy.assumptions.assume import global_assumptions
import sympy as sym
import re
from sympy.utilities.lambdify import lambdify
import sys
from copy import copy


# Note imports at bottom to avoid circular dependencies

# TODO, propagate assumptions for arithmetic.........  This may be
# tricky.  At the moment only a limited propagation of assumptions are
# performed.

__all__ = ('pprint', 'pretty', 'latex', 'DeltaWye', 'WyeDelta', 'tf',
           'symbol', 'sympify',
           'zp2tf', 'Expr', 's', 'sExpr', 't', 'tExpr', 'f', 'fExpr', 'cExpr',
           'omega', 'omega1', 'omegaExpr', 'Phasor',
           'pi', 'cos', 'sin', 'tan', 'atan', 'atan2',
           'exp', 'sqrt', 'log', 'log10', 'gcd', 'oo', 'inf',
           'H', 'Heaviside', 'DiracDelta', 'j', 'u', 'delta',
           'Vector', 'Matrix', 'VsVector', 'IsVector', 'YsVector', 'ZsVector',
           'Hs', 'Is', 'Vs', 'Ys', 'Zs',
           'Ht', 'It', 'Vt', 'Yt', 'Zt',
           'Hf', 'If', 'Vf', 'Yf', 'Zf',
           'Iphasor', 'Vphasor', 'Yphasor', 'Zphasor',
           'Homega', 'Iomega', 'Vomega', 'Yomega', 'Zomega')

func_pattern = re.compile(r"\\operatorname{(.*)}")

all_assumptions = ('ac', 'dc', 'causal')

from sympy.printing.latex import LatexPrinter 
from sympy.printing.pretty.pretty import PrettyPrinter 

init = False

class LcapyLatexPrinter(LatexPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr

        if expr == sym.I:
            return "j"
        return super(LcapyLatexPrinter, self)._print(expr)


def latex(expr, **settings):
    return LcapyLatexPrinter(settings).doprint(expr)


class LcapyPrettyPrinter(PrettyPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr

        if expr == sym.I:
            return self._print_basestring("j")
        return super(LcapyPrettyPrinter, self)._print(expr)


def pretty(expr, **settings):
    return LcapyPrettyPrinter(settings).doprint(expr)


def uppercase_name(name):

    return name[0].upper() + name[1:]


class Context(object):

    def __init__(self):
        self.symbols = {}
        self.assumptions = {}
        self.previous = None

    def new(self):

        new_context = Context()
        new_context.symbols.update(self.symbols)
        new_context.assumptions.update(self.assumptions)
        return new_context

    def switch(self):

        global context

        self.previous = context
        context = self
        global_assumptions.clear()
        global_assumptions.update(self.assumptions)

    def restore(self):
        
        if self.previous is None:
            return

        self.assumptions.update(global_assumptions)
        global_assumptions.clear()
        global_assumptions.update(self.previous.assumptions)


def sympify(expr, evaluate=True, **assumptions):
    """Create a sympy expression."""

    return sympify1(expr, context.symbols, 
                    evaluate, **assumptions)

def symbol(name, **assumptions):

    return sympify(name, **assumptions)


global_context = Context()
context = global_context

ssym = symbol('s')
tsym = symbol('t', real=True)
fsym = symbol('f', real=True)
omegasym = symbol('omega', real=True)
omega1sym = symbol('omega_1', real=True)


class Exprdict(dict):

    """Decorator class for dictionary created by sympy"""

    def pprint(self):
        """Pretty print"""

        return pprint(self)

    def latex(self):
        """Latex"""

        return latex_str(latex(self))


class Expr(object):

    """Decorator class for sympy classes derived from sympy.Expr"""

    # Perhaps have lookup table for operands to determine
    # the resultant type?  For example, Vs / Vs -> Hs
    # Vs / Is -> Zs,  Is * Zs -> Vs

    def __init__(self, arg, **assumptions):

        if isinstance(arg, Expr):
            if assumptions == {}:
                assumptions = arg.assumptions
            arg = arg.expr

        # Perhaps could set dc.
        if arg == 0:
            assumptions['causal'] = True

        # There are two types of assumptions.
        #   1. There are the sympy assumptions that are only associated
        #      with symbols, for example, real=True.
        #   2. The expr assumptions such as dc, ac, causal.

        for attr in all_assumptions:
            setattr(self, 'is_' + attr, assumptions.pop(attr, None))
        self.expr = sympify(arg, **assumptions)

    @property
    def val(self):
        """Return floating point value of expression if it can be evaluated,
        otherwise the expression"""

        return self.evalf()

    def __hash__(self):
        # This is needed for Python3 so can create a dict key,
        # say for subs.
        return hash(self.expr)

    def __getattr__(self, attr):

        # This gets called if there is no explicit attribute attr for
        # this instance.  We call the method of the wrapped sympy
        # class and rewrap the returned value if it is a sympy Expr
        # object.

        expr = self.expr
        if hasattr(expr, attr):
            a = getattr(expr, attr)

            # If it is not callable, directly wrap it.
            if not hasattr(a, '__call__'):
                if not isinstance(a, sym.Expr):
                    return a
                return self.__class__(ret)

            # If it is callable, create a function to pass arguments
            # through and wrap its return value.
            def wrap(*args):
                ret = a(*args)

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

        p.text(pretty(self.expr))

    def _repr_latex_(self):

        return '$%s$' % latex_str(self.latex())

    @property
    def assumptions(self):

        assumptions = {}
        for attr in all_assumptions:
            assumption = getattr(self, 'is_' + attr)
            if assumption is not None:
                assumptions[attr] = assumption
        return assumptions

    @assumptions.setter
    def assumptions(self, **assumptions):

        for attr in all_assumptions:
            if attr in assumptions:
                setattr(self, 'is_' + attr, assumptions.pop(attr))
        if assumptions != {}:
            raise ValueError('Unknown assumption %s' % assumptions)

    def __abs__(self):
        """Absolute value"""

        return self.__class__(abs(self.expr), **self.assumptions)

    def __neg__(self):
        """Negation"""

        return self.__class__(-self.expr, **self.assumptions)

    def __compat__(self, x):
        """Check if args are compatible and if so return compatible class."""

        # Could disallow Vs + Is, etc.
        # Could also convert Vs / Zs -> Is, etc.
        # But, what about (Vs * Vs) / (Vs * Is) ???

        cls = self.__class__
        if not isinstance(x, Expr):
            return cls

        xcls = x.__class__
        if cls == xcls:
            return cls

        if xcls in (Expr, cExpr):
            return cls

        if cls in (Expr, cExpr):
            return xcls

        if isinstance(x, cls):
            return xcls

        if isinstance(self, xcls):
            return cls

        if isinstance(self, tExpr) and isinstance(x, tExpr):
            return cls

        if isinstance(self, sExpr) and isinstance(x, sExpr):
            return cls

        if isinstance(self, Phasor) and isinstance(x, Phasor):
            return cls

        if isinstance(self, Phasor) and isinstance(x, Expr):
            return cls
        
        raise ValueError('Cannot combine %s(%s) with %s(%s)' % 
                         (cls.__name__, self, xcls.__name__, x))


    def __compat_add__(self, x):

        cls = self.__class__
        if not isinstance(x, Expr):
            return cls, self, cls(x)

        xcls = x.__class__
        if cls == xcls:
            return cls, self, x

        if isinstance(self, Phasor) and isinstance(x, Expr):
            self = self.laplace()
            cls = self.__class__
            return cls, self, cls(x)

        if isinstance(self, Expr) and isinstance(x, Phasor):
            x = x.laplace()
            cls = self.__class__
            return cls, self, x

        cls = self.__compat__(x)
        return cls, self, cls(x)

    def __rdiv__(self, x):
        """Reverse divide"""

        cls = self.__compat__(x)
        x = cls(x)
        return cls(x.expr / self.expr)

    def __rtruediv__(self, x):
        """Reverse true divide"""

        cls = self.__compat__(x)
        x = cls(x)
        return cls(x.expr / self.expr)

    def __mul__(self, x):
        """Multiply"""

        cls = self.__compat__(x)
        x = cls(x)
        return cls(self.expr * x.expr)

    def __rmul__(self, x):
        """Reverse multiply"""

        cls = self.__compat__(x)
        x = cls(x)
        return cls(self.expr * x.expr)

    def __div__(self, x):
        """Divide"""

        cls = self.__compat__(x)
        x = cls(x)
        return cls(self.expr / x.expr)

    def __truediv__(self, x):
        """True divide"""

        cls = self.__compat__(x)
        x = cls(x)
        return cls(self.expr / x.expr)

    def __add__(self, x):
        """Add"""

        cls, self, x = self.__compat_add__(x)
        return cls(self.expr + x.expr)

    def __radd__(self, x):
        """Reverse add"""

        cls, self, x = self.__compat_add__(x)
        return cls(self.expr + x.expr)

    def __rsub__(self, x):
        """Reverse subtract"""

        cls, self, x = self.__compat_add__(x)
        return cls(x.expr - self.expr)

    def __sub__(self, x):
        """Subtract"""

        cls, self, x = self.__compat_add__(x)
        return cls(self.expr - x.expr)

    def __pow__(self, x):
        """Pow"""

        cls = self.__compat__(x)
        x = cls(x)
        return cls(self.expr ** x.expr)

    def __or__(self, x):
        """Parallel combination"""

        return self.parallel(x)

    def __eq__(self, x):
        """Equality"""

        if x is None:
            return False

        cls = self.__compat__(x)
        x = cls(x)

        # This fails if one of the operands has the is_real attribute
        # end the other doesn't...
        return self.expr == x.expr

    def __ne__(self, x):
        """Inequality"""

        if x is None:
            return True

        cls = self.__compat__(x)
        x = cls(x)

        return self.expr != x.expr

    def parallel(self, x):
        """Parallel combination"""

        cls = self.__compat__(x)
        x = cls(x)
        return cls(self.expr * x.expr / (self.expr + x.expr))

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

        # Perhaps replace _{xxx} with _{\mathrm{xxx}} if len(xxx) > 1
        # We catch this later on with latex_str

        return string

    def pretty(self):
        """Make pretty string"""
        return pretty(self.expr)

    def prettyans(self, name):
        """Make pretty string with LHS name"""

        return pretty(sym.Eq(sympify(name), self.expr))

    def pprint(self):
        """Pretty print"""

        # If interactive use pretty, otherwise use latex
        if hasattr(sys, 'ps1'):
            print(self.pretty())
        else:
            print(latex_str(self.latex()))

    def pprintans(self, name):
        """Pretty print string with LHS name"""
        print(self.prettyans(name))

    def latex(self):
        """Make latex string"""

        string = latex(self.expr)
        match = func_pattern.match(string)
        if match is not None:
            # v_1(t) -> \operatorname{v\_1}\left( t \right)
            # operatorname requires amsmath so switch to mathrm
            string = r'\mathrm{%s}' % match.groups()[0].replace('\\_', '_')

        # sympy uses theta for Heaviside
        string = string.replace(r'\theta\left', r'u\left')
        return latex_str(string)

    def latexans(self, name):
        """Print latex string with LHS name"""

        expr = sym.Eq(sympify(name), self.expr)

        return latex_str(latex(expr))

    @property
    def N(self):
        """Return numerator of rational function"""

        return self.numerator

    @property
    def numerator(self):
        """Return numerator of rational function"""

        expr = self.expr
        if not expr.is_rational_function(self):
            raise ValueError('Expression not a rational function')

        numer, denom = expr.as_numer_denom()
        return self.__class__(numer)

    @property
    def D(self):
        """Return denominator of rational function"""

        return self.denominator

    @property
    def denominator(self):
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

        dst = self.__class__(sym.re(self.expr).simplify())
        dst.part = 'real'
        return dst

    @property
    def imag(self):
        """Return imaginary part"""

        dst = self.__class__(sym.im(self.expr).simplify())
        dst.part = 'imaginary'
        return dst

    @property
    def real_imag(self):
        """Rewrite as x + j * y"""
        
        return self.real + j * self.imag
    
    def rationalize_denominator(self):
        """Rationalize denominator by multiplying numerator and denominator by
        complex conjugate of denominator"""

        N = self.N
        D = self.D
        Dconj = D.conjugate
        Nnew = (N * Dconj).simplify()
        #Dnew = (D * Dconj).simplify()
        Dnew = (D.real**2 + D.imag**2).simplify()

        Nnew = Nnew.real_imag

        return Nnew / Dnew

    @property
    def magnitude(self):
        """Return magnitude"""

        R = self.rationalize_denominator()
        N = R.N
        Dnew = R.D
        Nnew = sqrt((N.real**2 + N.imag**2).simplify())

        dst = Nnew / Dnew
        dst.part = 'magnitude'
        return dst

    @property
    def abs(self):
        """Return magnitude"""

        return self.magnitude

    @property
    def dB(self):
        """Return magnitude in dB"""

        # Need to clip for a desired dynamic range?
        # Assume reference is 1. 
        dst = 20 * log10(self.magnitude)
        dst.part = 'magnitude'
        dst.units = 'dB'
        return dst

    @property
    def phase(self):
        """Return phase in radians"""

        R = self.rationalize_denominator()
        N = R.N
        
        G = gcd(N.real, N.imag)
        new = N / G
        dst = atan2(new.imag, new.real)
        dst.part = 'phase'
        dst.units = 'rad'
        return dst

    @property
    def phase_degrees(self):
        """Return phase in degrees"""

        dst = self.phase * 180.0 / pi
        dst.part = 'phase'
        dst.units = 'degrees'
        return dst

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

    def evaluate(self, arg=None):
        """Evaluate expression at arg.  arg may be a scalar, or a vector.

        There can be no symbols in the expression except for the variable.

        Note, expressions such as exp(-alpha*t) * Heaviside(t) will
        not evaluate correctly since the exp will overflow for -t and
        produce an Inf.  When this is multiplied by 0 from the
        Heaviside function we get Nan. """

        # Use symbol names to avoid problems with symbols of the same
        # name with different assumptions.
        varname = self.var.name
        free_symbols = set([symbol.name for symbol in self.expr.free_symbols])
        if varname in free_symbols:
            free_symbols -= set((varname, ))
        if free_symbols != set():
            raise ValueError('Undefined symbols %s in expression %s' % (tuple(free_symbols), self.var))

        if arg is None:
            if self.expr.find(self.var) != set():
                raise ValueError('Need value to evaluate expression at')
            # The arg is irrelevant since the expression is a constant.
            arg = 0

        # Perhaps should check if expr.args[1] == Heaviside('t') and not
        # evaluate if t < 0?

        def dirac(arg):

            return np.inf if arg == 0.0 else 0.0

        def sqrt(arg):

            if arg < 0:
                return 1j * np.sqrt(-arg)
            return np.sqrt(arg)

        # For negative arguments, np.sqrt will return Nan.
        # np.lib.scimath.sqrt converts to complex but cannot be used
        # for lamdification!
        func = lambdify(self.var, self.expr,
                        ({'DiracDelta' : dirac,
                          'sqrt' : sqrt},
                         "numpy", "sympy", "math"))

        if np.isscalar(arg):
            v1 = arg
        else:
            v1 = arg[0]

        try:
            result = func(v1)
            response = complex(result)
        except NameError:
            raise RuntimeError('Cannot evaluate expression %s' % self)
        except (AttributeError, TypeError):
            if self.expr.is_Piecewise:
                raise RuntimeError(
                    'Cannot evaluate expression %s,'
                    ' due to undetermined conditional result' % self)

            raise RuntimeError(
                'Cannot evaluate expression %s,'
                ' probably have a mysterious function' % self)

        if np.isscalar(arg):
            if np.allclose(response.imag, 0.0):
                response = response.real
            return response

        try:
            response = np.array([complex(func(v1)) for v1 in arg])
        except TypeError:
            raise TypeError(
                'Cannot evaluate expression %s,'
                ' probably have undefined symbols' % self)

        if np.allclose(response.imag, 0.0):
            response = response.real
        return response

    def _subs1(self, old, new, **kwargs):

        # This will fail is a variable has different attributes,
        # such as positive or real.
        # Should check for bogus substitutions, such as t for s.

        expr = new
        if isinstance(new, Expr):
            cls = new.__class__
            expr = new.expr
        else:
            cls = self.__class__
            expr = sympify(expr)

        class_map = {(Hs, omegaExpr) : Homega,
                     (Is, omegaExpr) : Iomega,
                     (Vs, omegaExpr) : Vomega,
                     (Ys, omegaExpr) : Yomega,
                     (Zs, omegaExpr) : Zomega,
                     (Hs, fExpr) : Hf,
                     (Is, fExpr) : If,
                     (Vs, fExpr) : Vf,
                     (Ys, fExpr) : Yf,
                     (Zs, fExpr) : Zf,
                     (Hf, omegaExpr) : Homega,
                     (If, omegaExpr) : Iomega,
                     (Vf, omegaExpr) : Vomega,
                     (Yf, omegaExpr) : Yomega,
                     (Zf, omegaExpr) : Zomega,
                     (Homega, fExpr) : Hf,
                     (Iomega, fExpr) : If,
                     (Vomega, fExpr) : Vf,
                     (Yomega, fExpr) : Yf,
                     (Zomega, fExpr) : Zf}

        if (self.__class__, new.__class__) in class_map:
            cls = class_map[(self.__class__, new.__class__)]

        name = canonical_name(old)

        if not isinstance(name, str):
            name = str(name)

        # Replace symbol names with symbol definitions to
        # avoid problems with real or positive attributes.
        if name not in context.symbols:
            raise ValueError('Unknown symbol %s' % old)
        old = context.symbols[name]

        return cls(self.expr.subs(old, expr))

    def __call__(self, arg):
        """Substitute arg for variable.  If arg is an tuple or list
        return a list.  If arg is an numpy array, return
        numpy array.

        See also evaluate.
        """

        if isinstance(arg, (tuple, list)):
            return [self._subs1(self.var, arg1) for arg1 in arg]

        if isinstance(arg, np.ndarray):
            return np.array([self._subs1(self.var, arg1) for arg1 in arg])

        return self._subs1(self.var, arg)

    def subs(self, *args, **kwargs):
        """Substitute variables in expression, see sympy.subs for usage"""

        if len(args) > 2:
            raise ValueError('Too many arguments')
        if len(args) == 0:
            raise ValueError('No arguments')

        if len(args) == 2:
            return self._subs1(args[0], args[1])

        if  isinstance(args[0], dict):
            dst = self
            for key, val in args[0].items():
                dst = dst._subs1(key, val, **kwargs)

            return dst

        return self._subs1(self.var, args[0])

    @property
    def label(self):

        label = ''
        if hasattr(self, 'quantity'):
            label += self.quantity
            if hasattr(self, 'part'):
                label += ' ' + self.part
        else:
            if hasattr(self, 'part'):
                label += uppercase_name(self.part)
        if hasattr(self, 'units') and self.units != '':
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

    def differentiate(self, arg=None):

        if arg is None:
            arg = self.var
        return self.__class__(sym.diff(self.expr, arg))

    def diff(self, arg=None):

        return self.differentiate(arg)

    def debug(self):

        expr = self.expr
        print(expr)
        for symbol in expr.free_symbols:
            print('  %s: %s' % (symbol, symbol.assumptions0))

class sfwExpr(Expr):

    def __init__(self, val, **assumptions):

        super(sfwExpr, self).__init__(val, **assumptions)

    def _as_ratfun_delay(self):
        """Split expr as (N, D, delay)
        where expr = (N / D) * exp(var * delay)
        
        Note, delay only represents a delay when var is s."""

        expr, var = self.expr, self.var

        delay = sympify(0)

        if expr.is_rational_function():
            numer, denom = expr.as_numer_denom()
            N = sym.Poly(numer, var)
            D = sym.Poly(denom, var)

            return N, D, delay

        # Note, there is a bug in sympy factor, TODO warn if detected.
        F = sym.factor(expr).as_ordered_factors()

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

        try:
            N, D, delay = self._as_ratfun_delay()
        except ValueError:
            return self.__class__(self.expr, **self.assumptions)

        K = sym.cancel(N.LC() / D.LC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        # Divide by leading coefficient
        N = N.monic()
        D = D.monic()

        expr = K * (N / D)

        return self.__class__(expr, **self.assumptions)

    def general(self):
        """Convert rational function to general form.

        See also canonical, partfrac, mixedfrac, and ZPK"""

        N, D, delay = self._as_ratfun_delay()

        expr = sym.cancel(N / D, self.var)
        if delay != 0:
            expr *= sym.exp(self.var * delay)

        return self.__class__(expr, **self.assumptions)

    def partfrac(self):
        """Convert rational function into partial fraction form.

        See also canonical, mixedfrac, general, and ZPK"""

        F, R, Q, delay = self._as_residue_parts()

        expr = Q.as_expr()
        for f, r in zip(F, R):
            expr += r / f

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return self.__class__(expr, **self.assumptions)

    def mixedfrac(self):
        """Convert rational function into mixed fraction form.

        See also canonical, general, partfrac and ZPK"""

        N, D, delay = self._as_ratfun_delay()

        Q, M = N.div(D)

        expr = Q + M / D

        if delay != 0:
            expr *= sym.exp(-self.var * delay)

        return self.__class__(expr, **self.assumptions)

    def ZPK(self):
        """Convert to pole-zero-gain (PZK) form.

        See also canonical, general, mixedfrac, and partfrac"""

        N, D, delay = self._as_ratfun_delay()

        K = sym.cancel(N.LC() / D.LC())
        if delay != 0:
            K *= sym.exp(self.var * delay)

        zeros = sym.roots(N)
        poles = sym.roots(D)

        return self.__class__(_zp2tf(zeros, poles, K, self.var), **self.assumptions)


class sExpr(sfwExpr):

    """s-domain expression or symbol"""

    var = ssym

    def __init__(self, val, **assumptions):

        super(sExpr, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = tExpr

        if self.expr.find(tsym) != set():
            raise ValueError(
                's-domain expression %s cannot depend on t' % self.expr)

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
        return self(sym.I * w)

    def initial_value(self):
        """Determine value at t = 0"""

        return self.__class__(sym.limit(self.expr * self.var, self.var, sym.oo))

    def final_value(self):
        """Determine value at t = oo"""

        return self.__class__(sym.limit(self.expr * self.var, self.var, 0))

    def _inverse_laplace(self, **assumptions):

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
        for factor in expr.as_ordered_factors():
            if factor == sym.oo:
                return factor

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
                    # conjugate.  Unforrtunately, for symbolic expressions
                    # we cannot tell if a quadratic has two real poles,
                    # a repeat real pole, or a complex conjugate pair of poles.
                    P2[pc] = 0
                  
                    p_re = sym.re(p)
                    p_im = sym.im(p)
                    r_re = sym.re(r)
                    r_im = sym.im(r)
                    etd = sym.exp(p_re * td)
                    result2 += 2 * r_re * etd * sym.cos(p_im * td)
                    result2 -= 2 * r_im * etd * sym.sin(p_im * td)
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

        if assumptions.get('ac', False):
            return result1 + result2

        if assumptions.get('causal', False):
            return result1 + result2 * sym.Heaviside(td)

        if delay != 0:
            result2 *= sym.Heaviside(td)
                
        return sym.Piecewise((result1 + result2, 't>=0'))


    def inverse_laplace(self, **assumptions):
        """Attempt inverse Laplace transform.

        If causal=True the response is zero for t < 0 and
        the result is multiplied by Heaviside(t)
        If ac=True or dc=True the result is extrapolated for t < 0.
        Otherwise the result is only known for t >= 0.

        """

        if assumptions == {}:
            assumptions = self.assumptions

        if assumptions.get('dc', False):
            result = self * s
            
            free_symbols = set([symbol.name for symbol in result.free_symbols])
            if 's' in free_symbols:
                raise ValueError('Something wonky going on, expecting dc.'
                                 ' Perhaps have capacitors in series?')
            return self._laplace_conjugate_class(result)

        try:
            result = self._inverse_laplace(**assumptions)

        except:

            print('Determining inverse Laplace transform with sympy...')

            # Try splitting into partial fractions to help sympy.
            expr = self.partfrac().expr

            # This barfs when needing to generate Dirac deltas
            from sympy.integrals.transforms import inverse_laplace_transform
            result = inverse_laplace_transform(expr, t, self.var)

            if not assumptions.causal:
                result = sym.Piecewise((result, 't>=0'))

        if hasattr(self, '_laplace_conjugate_class'):
            result = self._laplace_conjugate_class(result)
        return result

    def time(self, **assumptions):
        return self.inverse_laplace(**assumptions)

    def laplace(self):
        """Convert to s-domain representation"""

        return self.__class__(self, **self.assumptions)

    def phasor(self, **assumptions):

        return self.time(**assumptions).phasor(**assumptions)

    def transient_response(self, tvector=None):
        """Evaluate transient (impulse) response"""

        texpr = self.inverse_laplace()

        if tvector is None:
            return texpr

        print('Evaluating inverse Laplace transform...')
        return texpr.evaluate(tvector)

    def impulse_response(self, tvector=None):
        """Evaluate transient (impulse) response"""

        return self.transient_response(tvector)

    def step_response(self, tvector=None):
        """Evaluate step response"""

        return (self / self.var).transient_response(tvector)

    def angular_frequency_response(self, wvector=None):
        """Convert to angular frequency domain and evaluate response if
        angular frequency vector specified

        """

        X = self(j * omega)

        if wvector is None:
            return X

        return X.evaluate(wvector)

    def frequency_response(self, fvector=None):
        """Convert to frequency domain and evaluate response if frequency
        vector specified

        """

        X = self(j * 2 * pi * f)

        if fvector is None:
            return X

        return X.evaluate(fvector)

    def response(self, x, t):
        """Evaluate response to input signal x at times t"""

        if len(x) != len(t):
            raise ValueError('x must have same length as t')

        dt = t[1] - t[0]
        if not np.allclose(np.diff(t), np.ones(len(t) - 1) * dt):
            raise (ValueError, 't values not equally spaced')

        N, D, delay = self._as_ratfun_delay()

        Q, M = N.div(D)
        expr = M / D

        N = len(t)

        # Evaluate transient response.
        th = np.arange(N) * dt - dt
        h = sExpr(expr).transient_response(th)

        print('Convolving...')
        ty = t
        y = np.convolve(x, h)[0:N] * dt

        if Q:
            # Handle Dirac deltas and their derivatives.
            C = Q.all_coeffs()
            for n, c in enumerate(C):

                y += c * x

                x = np.diff(x) / dt
                x = np.hstack((x, 0))

        from scipy.interpolate import interp1d

        if delay != 0.0:
            print('Interpolating...')
            # Try linear interpolation; should oversample first...
            y = interp1d(ty, y, bounds_error=False, fill_value=0)
            y = y(t - delay)

        return y

    def decompose(self):

        N, D, delay = self._as_ratfun_delay()

        return N, D, delay

    def evaluate(self, svector):

        return super(sExpr, self).evaluate(svector)

    def plot(self, t=None, **kwargs):
        """Plot pole-zero map"""

        from lcapy.plot import plot_pole_zero

        plot_pole_zero(self, **kwargs)


class tsExpr(sExpr):

    """t or s-domain expression or symbol, interpreted in time domain
    if not containing s"""

    def __init__(self, val):

        assumptions = {}

        # If no s in expression evaluate as tExpr and convert to s-domain.
        if 's' not in symbols_find(val):
            tval = tExpr(val)
            val = tval.laplace().expr
            val = val.subs(context.symbols)
            assumptions = tval.assumptions

        super(tsExpr, self).__init__(val, real=True, **assumptions)


class fExpr(sfwExpr):

    """Fourier domain expression or symbol"""

    var = fsym
    domain_name = 'Frequency'
    domain_units = 'Hz'

    def __init__(self, val):

        super(fExpr, self).__init__(val, real=True)
        # Define when class defined.
        self._fourier_conjugate_class = tExpr

        if self.expr.find(ssym) != set():
            raise ValueError(
                'f-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find(tsym) != set():
            raise ValueError(
                'f-domain expression %s cannot depend on t' % self.expr)

    def inverse_fourier(self):
        """Attempt inverse Fourier transform"""

        from sympy.integrals.transforms import inverse_fourier_transform

        result = inverse_fourier_transform(self.expr, tsym, self.var.expr)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        return result

    def plot(self, fvector=None, **kwargs):
        """Plot frequency response at values specified by fvector.
        
        There are many plotting options, see matplotlib.pyplot.plot.

        For example:
            V.plot(fvector, log_scale=True)
            V.real.plot(fvector, color='black')
            V.phase.plot(fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.        
        """

        from lcapy.plot import plot_frequency
        plot_frequency(self, fvector, **kwargs)


class omegaExpr(sfwExpr):

    """Fourier domain expression or symbol (angular frequency)"""

    var = omegasym
    domain_name = 'Angular frequency'
    domain_units = 'rad/s'

    def __init__(self, val):

        super(omegaExpr, self).__init__(val, real=True)
        self._fourier_conjugate_class = tExpr

        if self.expr.find(ssym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on s' % self.expr)
        if self.expr.find(tsym) != set():
            raise ValueError(
                'omega-domain expression %s cannot depend on t' % self.expr)

    def inverse_fourier(self):
        """Attempt inverse Fourier transform"""

        from sympy.integrals.transforms import inverse_fourier_transform

        result = inverse_fourier_transform(self.expr, tsym, self.var.expr)
        if hasattr(self, '_fourier_conjugate_class'):
            result = self._fourier_conjugate_class(result)
        return result

    def plot(self, wvector=None, **kwargs):
        """Plot angular frequency response at values specified by wvector.
        
        There are many plotting options, see matplotlib.pyplot.plot.

        For example:
            V.plot(fvector, log_scale=True)
            V.real.plot(fvector, color='black')
            V.phase.plot(fvector, color='black', linestyle='--')

        By default complex data is plotted as separate plots of magnitude (dB)
        and phase.        
        """

        from lcapy.plot import plot_angular_frequency
        plot_angular_frequency(self, wvector, **kwargs)


class tExpr(Expr):

    """t-domain expression or symbol"""

    var = tsym
    domain_name = 'Time'
    domain_units = 's'

    def __init__(self, val):

        super(tExpr, self).__init__(val, real=True)

        # Hack to avoid circular dep.
        if init:
            if is_dc(self, t):
                self.is_dc = True
            elif is_ac(self, t):
                self.is_ac = True
            elif is_causal(self, t):
                self.is_causal = True

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
            var = sym.Symbol('t')

        try:
            F, a, cond = sym.laplace_transform(self, var, ssym)
        except TypeError as e:
            expr = self.expr
            if not isinstance(expr, sym.function.AppliedUndef):
                raise TypeError(e)
            if (expr.nargs != 1) or (expr.args[0] != tsym):
                raise TypeError(e)
            # Convert v(t) to V(s), etc.
            if expr.func.__name__ == 'u':
                print('Warning, use Heaviside(t) for unit step u(t)')
            name = uppercase_name(expr.func.__name__) + '(s)'
            return sExpr(name)

        if hasattr(self, '_laplace_conjugate_class'):
            F = self._laplace_conjugate_class(F)
        return F

    def fourier(self):
        """Attempt Fourier transform"""

        F = sym.fourier_transform(self.expr, self.var, f)

        if hasattr(self, '_fourier_conjugate_class'):
            F = self._fourier_conjugate_class(F)
        return F

    def phasor(self, **assumptions):

        check = ACChecker(self, t)
        if not check.is_ac:
            raise ValueError('Do not know how to convert %s to phasor' % self)
        if check.omega != omega1sym:            
            raise ValueError('Angular frequency %s not omega1' % check.omega)
        phasor = Phasor(check.amp * exp(j * check.phase))
        return phasor

    def plot(self, t=None, **kwargs):

        from lcapy.plot import plot_time
        plot_time(self, t, **kwargs)

    def canonical(self):
        return self.__class__(self)


class cExpr(Expr):

    """Constant real expression or symbol.

    If the expression is known to be positive, use cExpr(expr, positive=True)
    """

    def __init__(self, val, **assumptions):

        symbols = symbols_find(val)
        if 's' in symbols:
            raise ValueError(
                'constant expression %s cannot depend on s' % val)
        if 't' in symbols:
            raise ValueError(
                'constant expression %s cannot depend on t' % val)

        assumptions['real'] = True
        super(cExpr, self).__init__(val, **assumptions)


class Phasor(sfwExpr):

    var = omega1sym

    def time(self, **assumptions):
        """Convert to time domain representation"""

        return self._laplace_conjugate_class(self.real.expr * cos(self.var * t) + self.imag.expr * sin(self.var * t))

    def laplace(self):
        """Convert to Laplace domain representation"""

        return self.time().laplace()

    def phasor(self):
        """Convert to phasor representation"""
        return self.__class__(self, **self.assumptions)


class Vphasor(Phasor):

    def __init__(self, val, **assumptions):

        super(Vphasor, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Vt

    def cpt(self):

        v = self
        if v.is_number:
            return Vac(v.expr)

        # Need a combination of components.
        return self


class Iphasor(Phasor):

    def __init__(self, val, **assumptions):

        super(Iphasor, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = It
    
    def cpt(self):

        i = self
        if i.is_number:
            return Iac(i.expr)

        # Need a combination of components.
        return self


class Zphasor(sfwExpr):
    """Phase impedance aka complex impedance"""

    def __init__(self, val, **assumptions):

        super(Zphasor, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Zt


class Yphasor(sfwExpr):
    """Phase admittance aka complex admittance"""

    def __init__(self, val, **assumptions):

        super(Yphasor, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Yt



s = sExpr('s')
t = tExpr('t')
f = fExpr('f')
omega = omegaExpr('omega')
omega1 = omegaExpr('omega_1')
pi = sym.pi
j = sym.I
oo = sym.oo
inf = sym.oo


def pprint(expr):

    print(pretty(expr))


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

        return pprint(self)

    def latex(self):

        return latex(self)

    def _reformat(self, method):
        """Helper method for reformatting expression"""

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

    # TODO. There is probably a cunning way to automatically handle
    # the following.

    def inv(self):

        return self.__class__(sym.Matrix(self).inv())

    def det(self):

        return sym.Matrix(self).det()


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
        
    return uMul(K, *(zz + pp))


def zp2tf(zeros, poles, K=1, var=None):
    """Create a transfer function from lists of zeros and poles,
    and from a constant gain"""

    return Hs(_zp2tf(zeros, poles, K, var))


# Perhaps use a factory to create the following classes?

class Zs(sExpr):

    """s-domain impedance value"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val, **assumptions):

        super(Zs, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Zt

    @classmethod
    def C(cls, Cval):

        return cls(1 / (s * Cval))

    @classmethod
    def G(cls, Gval):

        return cls(1 / Gval)

    @classmethod
    def L(cls, Lval):

        return cls(s * Lval)

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

    def __init__(self, val, **assumptions):

        super(Ys, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = Yt

    @classmethod
    def C(cls, Cval):

        return cls(s * Cval)

    @classmethod
    def G(cls, Gval):

        return cls(Gval)

    @classmethod
    def L(cls, Lval):

        return cls(1 / (s * Lval))

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

    def __init__(self, val, **assumptions):

        super(Vs, self).__init__(val, **assumptions)
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

    def __init__(self, val, **assumptions):

        super(Is, self).__init__(val, **assumptions)
        self._laplace_conjugate_class = It

    def cpt(self):

        i = self * s

        if i.is_number:
            return Idc(i.expr)

        # Need a combination of components.
        return self


class Hs(sExpr):

    """s-domain ratio"""

    quantity = 's-ratio'
    units = ''

    def __init__(self, val, **assumptions):

        super(Hs, self).__init__(val, **assumptions)
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

    """f-domain voltage (units V/Hz)"""

    quantity = 'Voltage spectrum'
    units = 'V/Hz'

    def __init__(self, val):

        super(Vf, self).__init__(val)
        self._fourier_conjugate_class = Vt


class If(fExpr):

    """f-domain current (units A/Hz)"""

    quantity = 'Current spectrum'
    units = 'A/Hz'

    def __init__(self, val):

        super(If, self).__init__(val)
        self._fourier_conjugate_class = It


class Hf(fExpr):

    """f-domain transfer function response"""

    quantity = 'Transfer function'
    units = ''

    def __init__(self, val):

        super(Hf, self).__init__(val)
        self._fourier_conjugate_class = Ht


class Yomega(omegaExpr):

    """omega-domain admittance"""

    quantity = 'Admittance'
    units = 'siemens'

    def __init__(self, val):

        super(Yomega, self).__init__(val)
        self._fourier_conjugate_class = Yt


class Zomega(omegaExpr):

    """omega-domain impedance"""

    quantity = 'Impedance'
    units = 'ohms'

    def __init__(self, val):

        super(Zomega, self).__init__(val)
        self._fourier_conjugate_class = Zt


class Vomega(omegaExpr):

    """omega-domain voltage (units V/rad/s)"""

    quantity = 'Voltage spectrum'
    units = 'V/rad/s'

    def __init__(self, val):

        super(Vomega, self).__init__(val)
        self._fourier_conjugate_class = Vt


class Iomega(omegaExpr):

    """omega-domain current (units A/rad/s)"""

    quantity = 'Current spectrum'
    units = 'A/rad/s'

    def __init__(self, val):

        super(Iomega, self).__init__(val)
        self._fourier_conjugate_class = It


class Homega(omegaExpr):

    """omega-domain transfer function response"""

    quantity = 'Transfer function'
    units = ''

    def __init__(self, val):

        super(Homega, self).__init__(val)
        self._fourier_conjugate_class = Ht


class VsVector(Vector):

    _typewrap = Vs


class IsVector(Vector):

    _typewrap = Is


class YsVector(Vector):

    _typewrap = Ys


class ZsVector(Vector):

    _typewrap = Zs


def _funcwrap(func, *args):

    cls = args[0].__class__

    tweak_args = list(args)
    for m, arg in enumerate(args):
        if isinstance(arg, Expr):
            tweak_args[m] = arg.expr

    return cls(func(*tweak_args))


def sin(expr):

    return _funcwrap(sym.sin, expr)


def cos(expr):

    return _funcwrap(sym.cos, expr)


def tan(expr):

    return _funcwrap(sym.tan, expr)


def atan(expr):

    return _funcwrap(sym.atan, expr)


def atan2(expr1, expr2):

    return _funcwrap(sym.atan2, expr1, expr2)


def gcd(expr1, expr2):

    return _funcwrap(sym.gcd, expr1, expr2)


def exp(expr):

    return _funcwrap(sym.exp, expr)


def sqrt(expr):

    return _funcwrap(sym.sqrt, expr)


def log(expr):

    return _funcwrap(sym.log, expr)


def log10(expr):

    return _funcwrap(sym.log, expr, 10)


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


def VV(val, **assumptions):

    if assumptions.get('ac', False):
        return Vphasor(val)
    else:
        return Vs(val, **assumptions).canonical()


def II(val, **assumptions):

    if assumptions.get('ac', False):
        return Iphasor(val)
    else:
        return Is(val, **assumptions).canonical()


init = True
from lcapy.oneport import L, C, R, G, Idc, Vdc, Iac, Vac
