"""This file provides the Expr class.  This attempts to create a
consistent interface to SymPy's expressions.

Copyright 2014--2019 Michael Hayes, UCECE

"""


# TODO, propagate assumptions for arithmetic.........  This may be
# tricky.  At the moment only a limited propagation of assumptions are
# performed.

from __future__ import division
from .latex import latex_str
from .acdc import is_dc, is_ac, is_causal
from .ratfun import Ratfun, _zp2tf
from .sym import sympify, symsimplify, j, omegasym
from .context import context
from .printing import pprint, pretty, print_str
from .sympify import canonical_name
from .functions import sqrt, log10, atan2, gcd
import numpy as np
import sympy as sym
from sympy.utilities.lambdify import lambdify
import sys
from copy import copy
import six


def uppercase_name(name):

    return name[0].upper() + name[1:]


class Exprdict(dict):

    """Decorator class for dictionary created by sympy."""

    def pprint(self):
        """Pretty print"""

        return pprint(self)

    def latex(self):
        """Latex"""

        return latex_str(latex(self))

    def _repr_pretty_(self, p, cycle):

        p.text(pretty(self))


class Expr(object):

    """Decorator class for sympy classes derived from sympy.Expr"""

    one_sided = False
    var = None

    # Perhaps have lookup table for operands to determine
    # the resultant type?  For example, Vs / Vs -> Hs
    # Vs / Is -> Zs,  Is * Zs -> Vs
    # But what about Vs**2 ?

    def __init__(self, arg, **assumptions):

        if isinstance(arg, Expr):
            if assumptions == {}:
                assumptions = arg.assumptions.copy()
            arg = arg.expr

        # Perhaps could set dc?
        if arg == 0:
            assumptions['causal'] = True

        # There are two types of assumptions.
        #   1. There are the sympy assumptions that are only associated
        #      with symbols, for example, real=True.
        #   2. The expr assumptions such as dc, ac, causal.  These
        #      are primarily to help the inverse Laplace transform for sExpr
        #      classes.  The omega assumption is required for Phasors.

        self.assumptions = assumptions.copy()
        assumptions.pop('nid', None)
        
        self.expr = sympify(arg, **assumptions)

    @property
    def causal(self):
        return self.is_causal
        
    @causal.setter
    def causal(self, value):
        self.assumptions['causal'] = value
        if value:
            self.assumptions['dc'] = False
            self.assumptions['ac'] = False
        
    def infer_assumptions(self):
        self.assumptions['dc'] = None
        self.assumptions['ac'] = None
        self.assumptions['causal'] = None

    @property
    def is_dc(self):
        if 'dc' not in self.assumptions:
            self.infer_assumptions()
        return self.assumptions['dc'] == True

    @property
    def is_ac(self):
        if 'ac' not in self.assumptions:
            self.infer_assumptions()
        return self.assumptions['ac'] == True

    @property
    def is_causal(self):
        if 'causal' not in self.assumptions:
            self.infer_assumptions()
        return self.assumptions['causal'] == True

    @property
    def is_complex(self):
        if 'complex' not in self.assumptions:
            return False
        return self.assumptions['complex'] == True

    @property
    def val(self):
        """Return floating point value of expression if it can be evaluated,
        otherwise the expression."""

        return self.evalf()

    @property
    def omega(self):
        """Return angular frequency."""

        if 'omega' not in self.assumptions:
            return omegasym
        return self.assumptions['omega']

    def __hash__(self):
        # This is needed for Python3 so can create a dict key,
        # say for subs.
        return hash(self.expr)

# This will allow sym.sympify to magically extract the sympy expression
# but it will also bypass our __rmul__, __radd__, etc. methods that get called
# when sympy punts.
#
#    def _sympy_(self):
#        # This is called from sym.sympify
#        return self.expr

    def __getattr__(self, attr):

        # This gets called if there is no explicit attribute attr for
        # this instance.  We call the method of the wrapped sympy
        # class and rewrap the returned value if it is a sympy Expr
        # object.

        # FIXME.  This propagates the assumptions.  There is a
        # possibility that the operation may violate them.

        expr = self.expr
        if hasattr(expr, attr):
            a = getattr(expr, attr)

            # If it is not callable, directly wrap it.
            if not hasattr(a, '__call__'):
                if not isinstance(a, sym.Expr):
                    return a
                ret = a(*args)                
                if hasattr(self, 'assumptions'):
                    return self.__class__(ret, **self.assumptions)
                return self.__class__(ret)

            # If it is callable, create a function to pass arguments
            # through and wrap its return value.
            def wrap(*args):
                """This is wrapper for a SymPy function.
                For help, see the SymPy documentation."""

                ret = a(*args)

                if not isinstance(ret, sym.Expr):
                    return ret

                # Wrap the return value
                if hasattr(self, 'assumptions'):
                    return self.__class__(ret, **self.assumptions)
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

        return print_str(self.expr)

    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.expr)

    def _repr_pretty_(self, p, cycle):

        p.text(pretty(self.expr))

    def _repr_latex_(self):

        return '$%s$' % latex_str(self.latex())

    def __abs__(self):
        """Absolute value."""

        return self.__class__(self.abs, **self.assumptions)

    def __neg__(self):
        """Negation."""

        return self.__class__(-self.expr, **self.assumptions)

    def __compat_mul__(self, x, op):
        """Check if args are compatible and if so return compatible class."""

        # Could also convert Vs / Zs -> Is, etc.
        # But, what about (Vs * Vs) / (Vs * Is) ???

        assumptions = {}
        
        cls = self.__class__
        if not isinstance(x, Expr):
            return cls, self, cls(x), assumptions

        xcls = x.__class__

        if isinstance(self, sExpr) and isinstance(x, sExpr):
            if self.is_causal or x.is_causal:
                assumptions = {'causal' : True}
            elif self.is_dc and x.is_dc:
                assumptions = self.assumptions
            elif self.is_ac and x.is_ac:
                assumptions = self.assumptions
            elif self.is_ac and x.is_dc:
                assumptions = {'ac' : True}
            elif self.is_dc and x.is_ac:
                assumptions = {'ac' : True}                

        if cls == xcls:
            return cls, self, cls(x), assumptions

        # Allow omega * t but treat as t expression.
        if isinstance(self, omegaExpr) and isinstance(x, tExpr):
            return xcls, self, x, assumptions
        if isinstance(self, tExpr) and isinstance(x, omegaExpr):
            return cls, self, x, assumptions                    
        
        if xcls in (Expr, cExpr):
            return cls, self, cls(x), assumptions

        if cls in (Expr, cExpr):
            return xcls, self, x, assumptions

        if isinstance(x, cls):
            return xcls, self, cls(x), assumptions

        if isinstance(self, xcls):
            return cls, self, cls(x), assumptions

        if isinstance(self, tExpr) and isinstance(x, tExpr):
            return cls, self, cls(x), assumptions

        if isinstance(self, sExpr) and isinstance(x, sExpr):
            return cls, self, cls(x), assumptions

        if isinstance(self, omegaExpr) and isinstance(x, omegaExpr):
            return cls, self, cls(x), assumptions

        raise ValueError('Cannot combine %s(%s) with %s(%s) for %s' %
                         (cls.__name__, self, xcls.__name__, x, op))

    def __compat_add__(self, x, op):

        # Disallow Vs + Is, etc.

        assumptions = {}

        cls = self.__class__
        if not isinstance(x, Expr):
            return cls, self, cls(x), assumptions

        xcls = x.__class__

        if isinstance(self, sExpr) and isinstance(x, sExpr):
            if self.assumptions == x.assumptions:
                assumptions = self.assumptions
        
        if cls == xcls:
            return cls, self, x, assumptions

        # Handle Vs + sExpr etc.
        if isinstance(self, xcls):
            return cls, self, x, assumptions

        # Handle sExpr + Vs etc.
        if isinstance(x, cls):
            return xcls, self, cls(x), assumptions

        if xcls in (Expr, cExpr):
            return cls, self, x, assumptions

        if cls in (Expr, cExpr):
            return xcls, cls(self), x, assumptions

        raise ValueError('Cannot combine %s(%s) with %s(%s) for %s' %
                         (cls.__name__, self, xcls.__name__, x, op))

    def __rdiv__(self, x):
        """Reverse divide"""

        cls, self, x, assumptions = self.__compat_mul__(x, '/')
        return cls(x.expr / self.expr, **assumptions)

    def __rtruediv__(self, x):
        """Reverse true divide"""

        cls, self, x, assumptions = self.__compat_mul__(x, '/')
        return cls(x.expr / self.expr, **assumptions)

    def __mul__(self, x):
        """Multiply"""
        from .sup import Super

        if isinstance(x, Super):
            return x.__mul__(self)
        
        cls, self, x, assumptions = self.__compat_mul__(x, '*')
        return cls(self.expr * x.expr, **assumptions)

    def __rmul__(self, x):
        """Reverse multiply"""

        cls, self, x, assumptions = self.__compat_mul__(x, '*')
        return cls(self.expr * x.expr, **assumptions)

    def __div__(self, x):
        """Divide"""

        cls, self, x, assumptions = self.__compat_mul__(x, '/')
        return cls(self.expr / x.expr, **assumptions)

    def __truediv__(self, x):
        """True divide"""

        cls, self, x, assumptions = self.__compat_mul__(x, '/')
        return cls(self.expr / x.expr, **assumptions)

    def __add__(self, x):
        """Add"""

        cls, self, x, assumptions = self.__compat_add__(x, '+')
        return cls(self.expr + x.expr, **assumptions)

    def __radd__(self, x):
        """Reverse add"""

        cls, self, x, assumptions = self.__compat_add__(x, '+')
        return cls(self.expr + x.expr, **assumptions)

    def __rsub__(self, x):
        """Reverse subtract"""

        cls, self, x, assumptions = self.__compat_add__(x, '-')
        return cls(x.expr - self.expr, **assumptions)

    def __sub__(self, x):
        """Subtract"""

        cls, self, x, assumptions = self.__compat_add__(x, '-')
        return cls(self.expr - x.expr, **assumptions)

    def __pow__(self, x):
        """Pow"""

        # TODO: FIXME
        cls, self, x, assumptions = self.__compat_mul__(x, '**')
        return cls(self.expr ** x.expr, **assumptions)

    def __or__(self, x):
        """Parallel combination"""

        return self.parallel(x)

    def __eq__(self, x):
        """Equality"""

        # Note, this is used by the in operator.

        if x is None:
            return False

        try:
            cls, self, x, assumptions = self.__compat_add__(x, '==')
        except ValueError:
            return False
            
        x = cls(x)

        # This fails if one of the operands has the is_real attribute
        # and the other doesn't...
        return self.expr == x.expr

    def __ne__(self, x):
        """Inequality"""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '!=')
        x = cls(x)

        return self.expr != x.expr

    def __gt__(self, x):
        """Greater than"""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '>')
        x = cls(x)

        return self.expr > x.expr

    def __ge__(self, x):
        """Greater than or equal"""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '>=')
        x = cls(x)

        return self.expr >= x.expr

    def __lt__(self, x):
        """Less than"""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '<')
        x = cls(x)

        return self.expr < x.expr

    def __le__(self, x):
        """Less than or equal"""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '<=')
        x = cls(x)

        return self.expr <= x.expr

    def parallel(self, x):
        """Parallel combination"""

        cls, self, x, assumptions = self.__compat_add__(x, '|')
        x = cls(x)

        return cls(self.expr * x.expr / (self.expr + x.expr), **assumptions)

    def _pretty(self, *args, **kwargs):
        """Make pretty string."""

        # This works in conjunction with Printer._print
        # It is a hack to allow printing of _Matrix types
        # and its elements.
        expr = self.expr
        printer = args[0]

        return printer._print(expr)

    def _latex(self, *args, **kwargs):
        """Make latex string."""

        # This works in conjunction with LatexPrinter._print
        # It is a hack to allow printing of _Matrix types
        # and its elements.
        expr = self.expr
        printer = args[0]

        string = printer._print(expr)
        # sympy uses theta for Heaviside , use u(t) although I prefer H(t)
        string = string.replace(r'\theta\left', r'u\left')

        return string

    def pretty(self):
        """Make pretty string."""
        return pretty(self.expr)

    def prettyans(self, name):
        """Make pretty string with LHS name."""

        return pretty(sym.Eq(sympify(name), self.expr))

    def pprint(self):
        """Pretty print"""
        pprint(self)

    def pprintans(self, name):
        """Pretty print string with LHS name."""
        print(self.prettyans(name))

    def latex(self):
        """Make latex string."""

        string = latex(self.expr)
        match = func_pattern.match(string)
        if match is not None:
            # v_1(t) -> \operatorname{v\_1}\left( t \right)
            # operatorname requires amsmath so switch to mathrm
            string = r'\mathrm{%s}' % match.groups()[0].replace('\\_', '_')

        return latex_str(string)

    def math_latex(self):
        """Make latex math-mode string."""

        return '$' + self.latex() + '$'

    def latexans(self, name):
        """Print latex string with LHS name."""

        expr = sym.Eq(sympify(name), self.expr)

        return latex_str(latex(expr))

    @property
    def conjugate(self):
        """Return complex conjugate."""

        return self.__class__(sym.conjugate(self.expr), **self.assumptions)

    @property
    def real(self):
        """Return real part."""

        assumptions = self.assumptions.copy()
        assumptions['real'] = True        

        dst = self.__class__(symsimplify(sym.re(self.expr)), **assumptions)
        dst.part = 'real'
        return dst

    @property
    def imag(self):
        """Return imaginary part."""

        assumptions = self.assumptions.copy()
        assumptions['real'] = True
        
        dst = self.__class__(symsimplify(sym.im(self.expr)), **assumptions)
        dst.part = 'imaginary'
        return dst

    @property
    def real_imag(self):
        """Rewrite as x + j * y"""

        return self.real + j * self.imag

    @property
    def _ratfun(self):
        if self.var is None:
            raise ValueError('Not a rational function')
        
        return Ratfun(self.expr, self.var)

    @property
    def N(self):
        """Return numerator of rational function."""

        return self.numerator

    @property
    def D(self):
        """Return denominator of rational function."""

        if self.var is None:
            return self.__class__(1)
        
        return self.denominator

    @property
    def numerator(self):
        """Return numerator of rational function."""

        if self.var is None:
            return self

        return self.__class__(self._ratfun.numerator)

    @property
    def denominator(self):
        """Return denominator of rational function."""

        return self.__class__(self._ratfun.denominator)

    def rationalize_denominator(self):
        """Rationalize denominator by multiplying numerator and denominator by
        complex conjugate of denominator."""

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

        if self.is_real:
            dst = self
        else:
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
    def sign(self):
        """Return sign"""

        return self.__class__(sym.sign(self.expr))

    @property
    def dB(self):
        """Return magnitude in dB."""

        # Need to clip for a desired dynamic range?
        # Assume reference is 1.
        dst = 20 * log10(self.magnitude)
        dst.part = 'magnitude'
        dst.units = 'dB'
        return dst

    @property
    def phase(self):
        """Return phase in radians."""

        R = self.rationalize_denominator()
        N = R.N

        if N.imag == 0:
            dst = N.imag
        else:
            if N.real != 0:
                G = gcd(N.real, N.imag)
                N = N / G
            dst = atan2(N.imag, N.real)
            
        dst.part = 'phase'
        dst.units = 'rad'
        return dst

    @property
    def phase_degrees(self):
        """Return phase in degrees."""

        dst = self.phase * 180.0 / sym.pi
        dst.part = 'phase'
        dst.units = 'degrees'
        return dst

    @property
    def angle(self):
        """Return phase angle"""

        return self.phase

    @property
    def is_number(self):

        return self.expr.is_number

    @property
    def is_constant(self):

        expr = self.expr

        # Workaround for sympy bug
        # a = sym.sympify('DiracDelta(t)')
        # a.is_constant()
        
        if expr.has(sym.DiracDelta):
            return False
        
        return expr.is_constant()

    def evaluate(self, arg=None):
        """Evaluate expression at arg.  arg may be a scalar, or a vector.
        The result is of type float or complex.

        There can be no symbols in the expression except for the variable.
        """

        def evaluate_expr(expr, var, arg):

            # For some reason the new lambdify will convert a float
            # argument to complex
            
            def exp(arg):

                # Hack to handle exp(-a * t) * Heaviside(t) for t < 0
                # by trying to avoid inf when number overflows float.

                if isinstance(arg, complex):
                    if arg.real > 500:
                        arg = 500 + 1j * arg.imag
                elif arg > 500:
                    arg = 500;                        

                return np.exp(arg)

            def dirac(arg):
                return np.inf if arg == 0.0 else 0.0

            def heaviside(arg):

                return 1.0 if arg >= 0.0 else 0.0

            def sqrt(arg):
                if not isinstance(arg, complex) and arg < 0:
                    arg = arg + 0j
                return np.sqrt(arg)

            try:
                arg0 = arg[0]
                scalar = False
            except:
                arg0 = arg
                scalar = True

            # For negative arguments, np.sqrt will return Nan.
            # np.lib.scimath.sqrt converts to complex but cannot be used
            # for lamdification!
            func = lambdify(var, expr,
                            ({'DiracDelta' : dirac,
                              'Heaviside' : heaviside,
                              'sqrt' : sqrt, 'exp' : exp},
                             "numpy", "sympy", "math"))

            try:
                result = func(arg0)
                response = complex(result)
            except NameError:
                raise RuntimeError('Cannot evaluate expression %s' % self)
            except (AttributeError, TypeError):
                if expr.is_Piecewise:
                    raise RuntimeError(
                        'Cannot evaluate expression %s,'
                        ' due to undetermined conditional result' % self)

                raise RuntimeError(
                    'Cannot evaluate expression %s,'
                    ' probably have a mysterious function' % self)

            if scalar:
                if np.allclose(response.imag, 0.0):
                    response = response.real
                return response

            try:
                response = np.array([complex(func(arg0)) for arg0 in arg])
            except TypeError:
                raise TypeError(
                    'Cannot evaluate expression %s,'
                    ' probably have undefined symbols' % self)

            if np.allclose(response.imag, 0.0):
                response = response.real
            return response

        expr = self.expr
        if hasattr(self, 'var'):
            var = self.var
            # Use symbol names to avoid problems with symbols of the same
            # name with different assumptions.
            varname = var.name
            free_symbols = set([symbol.name for symbol in expr.free_symbols])
            if varname in free_symbols:
                free_symbols -= set((varname, ))
            if free_symbols != set():
                raise ValueError('Undefined symbols %s in expression %s' % (tuple(free_symbols), self))

            if arg is None:
                if expr.find(var) != set():
                    raise ValueError('Need value to evaluate expression at')
                # The arg is irrelevant since the expression is a constant.
                arg = 0
        else:
            # Have no variable so must be a constant.
            var = None
            arg = 0

        try:
            arg = arg.evalf()
        except:
            pass

        if not (expr.is_Piecewise and expr.args[0].args[1] == (tsym >= 0)):            
            return evaluate_expr(expr, var, arg)

        try:
            arg0 = arg[0]
            scalar = False
        except:
            arg0 = arg
            scalar = True

        expr = expr.args[0].args[0]
            
        if scalar:
            if arg0 >= 0:
                return evaluate_expr(expr, var, arg)
            else:
                return sym.nan
        result =  evaluate_expr(expr, var, arg)
        result[arg < 0] = sym.nan
        return result

    def has(self, subexpr):
        """Test whether the sub-expression is contained."""
        if hasattr(subexpr, 'expr'):
            subexpr = subexpr.expr
        return self.expr.has(subexpr)

    def _subs1(self, old, new, **kwargs):

        # This will fail if a variable has different attributes,
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
        if name in context.symbols:
            old = context.symbols[name]
        else:
            # Perhaps have symbol defined using sympy?
            pass

        result = self.expr.subs(old, expr)

        # If get empty Piecewise, then result unknowable.  TODO: sympy
        # 1.2 requires Piecewise constructor to have at least one
        # pair.
        if False and result.is_Piecewise and result == sym.Piecewise():
            result = sym.nan

        # TODO: propagate assumptions?
        return cls(result)

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

    def limit(self, var, value, dir='+'):
        """Determine limit of expression(var) at var = value."""

        # Need to use lcapy sympify otherwise could use
        # getattr to call sym.limit.

        var = sympify(var)
        value = sympify(value)

        # Experimental.  Compare symbols by names.
        symbols = list(self.expr.free_symbols)
        symbolnames = [str(symbol) for symbol in symbols]
        if str(var) not in symbolnames:
            return self
        var = symbols[symbolnames.index(str(var))]
        
        ret = sym.limit(self.expr, var, value)
        if hasattr(self, 'assumptions'):
            return self.__class__(ret, **self.assumptions)
        return self.__class__(ret)

    def simplify(self):
        """Simplify expression."""
        
        ret = symsimplify(self.expr)
        return self.__class__(ret, **self.assumptions)

    def subs(self, *args, **kwargs):
        """Substitute variables in expression, see sympy.subs for usage."""

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

    def canonical(self):
        return self.__class__(self)


from .cexpr import Iconst, Vconst, cExpr        
from .fexpr import Hf, If, Vf, Yf, Zf, fExpr    
from .sexpr import Hs, Is, Vs, Ys, Zs, sExpr
from .texpr import Ht, It, Vt, Yt, Zt, tExpr
from .omegaexpr import Homega, Iomega, Vomega, Yomega, Zomega, omegaExpr
