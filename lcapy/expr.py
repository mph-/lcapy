# -*- coding: utf-8 -*-

"""This module provides the Expr class.  This attempts to create a
consistent interface to SymPy's expressions.

Copyright 2014--2021 Michael Hayes, UCECE

"""


# TODO, propagate assumptions for arithmetic.........  This may be
# tricky.  At the moment only a limited propagation of assumptions are
# performed.

from __future__ import division
from .assumptions import Assumptions
from .domains import UndefinedDomain
from .quantity import UndefinedQuantity
from .ratfun import Ratfun
from .sym import sympify, symsimplify, j, omegasym, symdebug, AppliedUndef
from .sym import capitalize_name, tsym, symsymbol, symbol_map
from .state import state
from .printing import pprint, pretty, print_str, latex
from .functions import sqrt, log10, atan2, gcd, exp, Function
from .units import units, u as uu, dB
from .utils import as_N_D, as_sum
from .state import state
import numpy as np
import sympy as sym
from sympy.utilities.lambdify import lambdify
from .sym import simplify
from .simplify import simplify_sin_cos
from collections import OrderedDict

__all__ = ('expr', 'symbol', 'symbols', 'deg', 'rad', 'degrees', 'radians')

class ExprPrint(object):

    @property
    def _pexpr(self):
        """Return expression for printing."""

        if not hasattr(self, 'expr'):
            return self            

        if state.show_units:
            if state.canonical_units:
                return self.expr_with_canonical_units
            else:
                return self.expr_with_units        
        else:
            return self.expr                
    
    def __repr__(self):
        """This is called by repr(expr).  It is used, e.g., when printing
        in the debugger."""
        
        return '%s(%s)' % (self.__class__.__name__, print_str(self._pexpr))

    def _repr_pretty_(self, p, cycle):
        """This is used by jupyter notebooks to display an expression using
        unicode.  It is also called by IPython when displaying an
        expression.""" 

        p.text(pretty(self._pexpr))

    # Note, _repr_latex_ is handled at the end of this file.
        
    def pretty(self, **kwargs):
        """Make pretty string."""
        return pretty(self._pexpr, **kwargs)

    def prettyans(self, name, **kwargs):
        """Make pretty string with LHS name."""
        return pretty(sym.Eq(sympify(name), self._pexpr), **kwargs)

    def pprint(self, **kwargs):
        """Pretty print"""
        pprint(self._pexpr, **kwargs)

    def pprintans(self, name, **kwargs):
        """Pretty print string with LHS name."""
        print(self.prettyans(name, **kwargs))

    def latex(self, **kwargs):
        """Make latex string."""
        return latex(self._pexpr, **kwargs)

    def latex_math(self, **kwargs):
        """Make latex math-mode string."""
        return '$' + self.latex(**kwargs) + '$'

    def latexans(self, name, **kwargs):
        """Print latex string with LHS name."""
        return latex(sym.Eq(sympify(name), self._pexpr), **kwargs)

    def srepr(self):
        return sym.repr(self)


class ExprContainer(object):    

    def evaluate(self):
        
        """Evaluate each element to convert to floating point."""        
        return self.__class__([v.evalf() for v in self])
    
    def simplify(self):
        """Simplify each element."""
        
        return self.__class__([simplify(v) for v in self])

    @property    
    def symbols(self):
        """Return dictionary of symbols in the expression keyed by name."""
        
        symbols = {}
        for expr in list(self):
            symbols.update(expr.symbols)
        return symbols

    
class ExprMisc(object):

    def pdb(self):
        """Enter the python debugger."""
        
        import pdb; pdb.set_trace()
        return self

        
class ExprDict(ExprPrint, ExprContainer, ExprMisc, OrderedDict):

    """Decorator class for dictionary created by sympy."""

    def evaluate(self):
        """Evaluate each element to convert to floating point.
        The keys are also converted if possible to handle
        dictionaries of poles/zeros."""

        new = self.__class__()
        for k, v in self.items():
            try:
                k = k.evalf()
            except:
                pass
            try:
                v = v.evalf()
            except:
                pass            
                
            new[k] = v
        return new

    def simplify(self):
        """Simplify each element but not the keys."""

        new = self.__class__()
        for k, v in self.items():
            new[k] = simplify(v)
        return new

    def subs(self, *args, **kwargs):
        """Substitute variables in expression, see sympy.subs for usage."""

        new = self.__class__()
        for k, v in self.items():
            try:
                k = k.subs(*args, **kwargs)
            except:
                pass
            try:
                v = v.subs(*args, **kwargs)
            except:
                pass            
                
            new[k] = v
        return new

    
class ExprList(ExprPrint, list, ExprContainer, ExprMisc):
    """Decorator class for list created by sympy."""

    # Have ExprPrint first so that its _repr__pretty_ is called
    # in preference to list's one.  Alternatively, add explicit
    # _repr_pretty_ method here.
    
    def __init__(self, iterable=None, evalf=False):

        if iterable is None:
            iterable = []
        
        eiterable = []
        for item in iterable:
            if evalf:
                try:
                    item = item.evalf()
                except:
                    pass
            else:
                item = expr(item)
            eiterable.append(item)
        
        super (ExprList, self).__init__(eiterable)

    def subs(self, *args, **kwargs):
        """Substitute variables in expression, see sympy.subs for usage."""
        
        return expr([e.subs(*args, **kwargs) for e in self])
        

class ExprTuple(ExprPrint, tuple, ExprContainer, ExprMisc):
    """Decorator class for tuple created by sympy."""

    # Tuples are immutable, need to use __new__
    def __new__(cls, iterable):

        eiterable = [expr(e) for e in iterable]
        return super (ExprTuple, cls).__new__(cls, eiterable)

    def subs(self, *args, **kwargs):
        """Substitute variables in expression, see sympy.subs for usage."""
        
        return expr((e.subs(*args, **kwargs) for e in self))

    
class Expr(UndefinedDomain, UndefinedQuantity, ExprPrint, ExprMisc):
    """Decorator class for sympy classes derived from sympy.Expr."""

    var = None

    _mul_mapping = {('voltage', 'admittance'): 'current',
                    ('current', 'impedance'): 'voltage',
                    ('voltage', 'transfer'): 'voltage',
                    ('current', 'transfer'): 'current',
                    ('transfer', 'transfer'): 'transfer',
                    ('voltage', 'constant'): 'voltage',
                    ('current', 'constant'): 'current',
                    ('admittance', 'constant'): 'admittance',
                    ('impedance', 'constant'): 'impedance',
                    ('transfer', 'constant'): 'transfer',
                    ('constant', 'constant'): 'constant',
                    ('voltage', 'voltage'): 'voltagesquared',
                    ('current', 'current'): 'currentsquared',
                    ('admittance', 'impedance'): 'transfer',                    
                    ('admittance', 'admittance'): 'admittancesquared',
                    ('impedance', 'impedance'): 'impedancesquared',
                    ('voltage', 'current'): 'power',
                    ('voltagesquared', 'admittance'): 'power',
                    ('currentsquared', 'impedance'): 'power',                    
                    ('impedancesquared', 'admittance'): 'impedance',
                    ('admittancesquared', 'impedance'): 'admittance',
                    ('power', 'impedance'): 'voltagesquared',
                    ('power', 'admittance'): 'currentsquared',
                    ('admittancesquared', 'constant'): 'admittancesquared',
                    ('impedancesquared', 'constant'): 'impedancesquared',
                    ('voltagesquared', 'constant'): 'voltagesquared',
                    ('currentsquared', 'constant'): 'currentsquared',
                    ('power', 'constant'): 'power'}
    
    _div_mapping = {('voltage', 'impedance'): 'current',
                    ('current', 'admittance'): 'voltage',
                    ('voltage', 'transfer'): 'voltage',
                    ('current', 'transfer'): 'current',
                    ('transfer', 'transfer'): 'transfer',
                    ('voltage', 'current'): 'impedance',
                    ('current', 'voltage'): 'admittance',                    
                    ('current', 'current'): 'transfer',
                    ('voltage', 'voltage'): 'transfer',
                    ('voltage', 'constant'): 'voltage',
                    ('current', 'constant'): 'current',
                    ('impedance', 'constant'): 'impedance',
                    ('admittance', 'constant'): 'admittance',
                    ('transfer', 'constant'): 'transfer',
                    ('constant', 'impedance'): 'admittance',
                    ('constant', 'admittance'): 'impedance',
                    ('constant', 'transfer'): 'transfer',
                    ('constant', 'constant'): 'constant',
                    ('impedance', 'impedance'): 'transfer',
                    ('admittance', 'admittance'): 'transfer',                                        
                    ('voltagesquared', 'voltage'): 'voltage',
                    ('currentsquared', 'current'): 'current',
                    ('admittancesquared', 'admittance'): 'admittance',
                    ('impedancesquared', 'impedance'): 'impedance',
                    ('power', 'current'): 'voltage',
                    ('power', 'voltage'): 'current',                    
                    ('power', 'admittance'): 'voltagesquared',
                    ('power', 'voltagesquared'): 'admittance',                    
                    ('power', 'impedance'): 'currentsquared',
                    ('power', 'currentsquared'): 'impedance',            
                    ('impedance', 'admittance'): 'impedancesquared',
                    ('admittance', 'impedance'): 'admittancesquared',
                    ('voltagesquared', 'impedance'): 'power',
                    ('currentsquared', 'admittance'): 'power',
                    ('voltagesquared', 'power'): 'impedance',
                    ('currentsquared', 'power'): 'admittance',
                    ('admittancesquared', 'constant'): 'admittancesquared',
                    ('impedancesquared', 'constant'): 'impedancesquared',
                    ('voltagesquared', 'constant'): 'voltagesquared',
                    ('currentsquared', 'constant'): 'currentsquared',
                    ('admittancesquared', 'admittancesquared'): 'transfer',
                    ('impedancesquared', 'impedancesquared'): 'transfer',
                    ('voltagesquared', 'voltagesquared'): 'transfer',
                    ('currentsquared', 'currentsquared'): 'transfer',
                    ('power', 'power'): 'transfer',
                    ('power', 'constant'): 'power'}
    
    # This needs to be larger than what sympy defines so
    # that the __rmul__, __radd__ methods get called.
    # Otherwise pi * t becomes a Mul rather than a TimeDomainExpression object.
    _op_priority = 1000
    

    @property
    def _pexpr(self):
        """Return expression for printing."""

        if not hasattr(self, 'expr'):
            return self
        
        if state.show_units:
            if state.canonical_units:
                return self.expr_with_canonical_units
            else:
                return self.expr_with_units                    
        else:
            return self.expr                

    def __init__(self, arg, **assumptions):
        """

         There are three types of assumptions:
           1. The sympy assumptions associated with symbols, for example,
              real=True.
           2. The expr assumptions such as dc, ac, causal.  These are primarily
              to help the inverse Laplace transform for LaplaceDomainExpression classes.
           3. Additional parameters such as nid for NoiseDomain and omega
              for PhasorDomain.

        """

        if isinstance(arg, Expr):
            ass = arg.assumptions.copy()
            if self.is_always_causal:
                ass.set('causal', True)
            
            self.assumptions = ass.merge(**assumptions)
            self.expr = arg.expr
            try:
                self._units = self._default_units
            except:
                self._units = sym.S.One
            return

        assumptions = Assumptions(assumptions)
        
        # Perhaps could set dc?
        #if arg == 0:
        #    assumptions.set('causal', True)
        if self.is_always_causal:
            assumptions.set('causal', True)        

        self.assumptions = assumptions
        # Remove Lcapy assumptions from SymPy expr.
        self.expr = sympify(arg, **self.assumptions.sympy_assumptions())
        try:
            self._units = self._default_units
        except:
            self._units = sym.S.One        


    def _class_by_quantity(self, quantity, domain=None):

        if domain is None:
            domain = self.domain
        return expressionclasses.get_quantity(domain, quantity)

    def _class_by_domain(self, domain):

        return expressionclasses.get_quantity(domain, self.quantity)    
        
    def as_quantity(self, quantity):

        if quantity == 'voltage':
            return self.as_voltage()
        elif quantity == 'current':
            return self.as_current()
        elif quantity == 'impedance':
            return self.as_impedance()
        elif quantity == 'admittance':
            return self.as_admittance()
        elif quantity == 'transfer':
            return self.as_transfer()
        elif quantity == 'power':
            return self.as_power()
        elif quantity == 'undefined':
            return self.as_expr()        
        raise ValueError('Unknown quantity %s for %s' % (quantity, self))

    def as_domain(self, domain):

        if domain == 'time':
            return self.as_time()
        elif domain == 'laplace':
            return self.as_laplace()
        elif domain == 'fourier':
            return self.as_fourier()
        elif domain == 'phasor':
            return self.as_phasor()        
        elif domain == 'angular fourier':
            return self.as_angular_fourier()
        raise ValueError('Unknown domain %s for %s' % (domain, self))

    def as_voltage(self):
        return self._class_by_quantity('voltage')(self)

    def as_current(self):
        return self._class_by_quantity('current')(self)

    def as_admittance(self):
        return self._class_by_quantity('admittance')(self)

    def as_impedance(self):
        return self._class_by_quantity('impedance')(self)

    def as_transfer(self):
        return self._class_by_quantity('transfer')(self)

    def as_power(self):
        return self._class_by_quantity('power')(self)

    def as_expr(self):
        return self

    def as_constant(self):
        if not self.is_unchanging:
            raise ValueError('Expression %s is not constant' % self)
        return self._class_by_quantity(self.quantity)(self)(cexpr(self))

    def as_superposition(self):
        from .superpositionvoltage import SuperpositionVoltage
        from .superpositioncurrent import SuperpositionCurrent
        
        if self.is_voltage:
            return SuperpositionVoltage(self)
        elif self.is_current:
            return SuperpositionCurrent(self)
        raise ValueError('Can only convert voltage or current to superposition')
    
    def as_time(self):
        return self.time()

    def as_laplace(self):
        return self.laplace()

    def as_phasor(self):
        return self.phasor()        

    def as_fourier(self):
        return self.fourier()

    def as_angular_fourier(self):
        return self.angular_fourier()
    
    def __str__(self, printer=None):
        """String representation of expression."""
        return print_str(self._pexpr)

    def __repr__(self):
        """This is called by repr(expr).  It is used, e.g., when printing
        in the debugger."""
        
        return '%s(%s)' % (self.__class__.__name__, print_str(self._pexpr))

    def _repr_pretty_(self, p, cycle):
        """This is used by jupyter notebooks to display an expression using
        unicode.  It is also called by IPython when displaying an
        expression.""" 

        p.text(pretty(self._pexpr))

    def _repr_latex_(self):
        """This is used by jupyter notebooks to display an expression using
        LaTeX markup.  However, this requires mathjax.  If this method
        is not defined, jupyter falls back on _repr_pretty_ which
        outputs unicode."""

        # This is called for Expr but not ExprList
        return '$$' + latex(self._pexpr) + '$$'        

    def _latex(self, *args, **kwargs):
        """Make latex string.  This is called by sympy.latex when it
        encounters an Expr type."""

        # This works in conjunction with LatexPrinter._print
        # It is a hack to allow printing of _Matrix types
        # and its elements.
        # This also catches sym.latex(expr) where expr is
        # an Lcapy expr.

        return self.latex(**kwargs)

    def _pretty(self, *args, **kwargs):
        """Make pretty string."""

        # This works in conjunction with Printer._print
        # It is a hack to allow printing of _Matrix types
        # and its elements.
        expr = self._pexpr
        printer = args[0]

        return printer._print(expr)

    @property
    def canonical_units(self):
        """Return the canonical units of the expression.  This is a simplified form,
        so volt * ampere becomes watt."""

        return units.simplify_units(self._units)        
    
    @property
    def units(self):
        """Return the units of the expression."""

        return self._units

    @units.setter
    def units(self, unit):
        """Set the units of the expression; these are simplified into canonical form."""        

        self._units = unit
    
    @property
    def is_causal(self):
        """Return True if zero for t < 0."""

        if self.assumptions.has_unspecified:
            self.assumptions.infer_from_expr(self)        
        return self.assumptions.is_causal

    @is_causal.setter
    def is_causal(self, value):

        self.assumptions.set('causal', value)

    @property
    def is_dc(self):

        if self.assumptions.has_unspecified:
            self.assumptions.infer_from_expr(self)                
        return self.assumptions.is_dc

    @is_dc.setter
    def is_dc(self, value):

        self.assumptions.set('dc', value)        

    @property
    def is_ac(self):

        if self.assumptions.has_unspecified:
            self.assumptions.infer_from_expr(self)                
        return self.assumptions.is_ac

    @is_ac.setter
    def is_ac(self, value):

        self.assumptions.set('ac', value)

    @property
    def is_unknown(self):
        """Return True if behaviour is unknown for t < 0."""

        if self.assumptions.has_unspecified:
            self.assumptions.infer_from_expr(self)                
        return self.assumptions.is_unknown

    @is_unknown.setter
    def is_unknown(self, value):

        self.assumptions.set('unknown', value)        

    @property
    def is_complex(self):
        
        if 'complex' not in self.assumptions:
            return False
        return self.assumptions['complex'] == True

    @property
    def is_conditional(self):
        """Return True if expression has a condition, such as t >= 0."""
        
        expr = self.expr
        # Could be more specific, such as self.var >= 0, but might
        # have self.var >= t1.
        return expr.is_Piecewise

    @property
    def is_rational_function(self):
        """Return True if expression is a rational function."""

        return self.expr.is_rational_function(self.var)
    
    @property
    def is_strictly_proper(self):
        """Return True if the degree of the dominator is greater
        than the degree of the numerator.
        This will throw an exception if the expression is not a
        rational function."""

        if self._ratfun is None:
            return False
        
        return self._ratfun.is_strictly_proper

    @property
    def is_phase(self):
        return self.part == 'phase'

    @property
    def is_phase_radians(self):
        return self.part == 'phase' and self.units == uu.rad

    @property
    def is_phase_degrees(self):
        return self.part == 'phase' and self.units == uu.deg

    @property
    def is_real_part(self):
        return self.part == 'real'

    @property
    def is_imag_part(self):
        return self.part == 'imaginary'

    @property
    def is_magnitude(self):
        return self.part == 'magnitude'

    @property
    def is_dB(self):
        return self.part == 'magnitude' and self.units == dB

    @property
    def ac(self):
        """Return the AC components."""

        return self.as_superposition().ac

    @property
    def dc(self):
        """Return the DC component."""

        return self.as_superposition().dc        

    @property
    def transient(self):
        """Return the transient component."""
        
        return self.as_superposition().transient
    
    @property
    def fval(self):
        """Evaluate expression and return as a python float value."""

        return float(self.val.expr)

    @property
    def cval(self):
        """Evaluate expression and return as a python complex value."""

        return complex(self.val.expr)
    
    @property
    def val(self):
        """Return floating point value of expression if it can be evaluated,
        otherwise the expression.

        This returns an Lcapy Expr object.   If you want a numerical value
        use expr.fval for a float value or expr.cval for a complex value."""

        return self.evalf()

    def evalf(self, n=15, *args, **kwargs):
        """Convert constants in an expression to floats, evaluated to `n`
        decimal places.  If the expression is a constant, return the
        floating point result.

        This returns an Lcapy Expr object.   If you want a numerical value
        use expr.fval for a float value or expr.cval for a complex value.

        See sympy.evalf for more details.

        """

        val = self.expr.evalf(n, *args, **kwargs)
        return self.__class__(val, **self.assumptions)    

    def __hash__(self):
        # This is needed for Python3 so can create a dict key,
        # say for subs.
        return hash(self.expr)


    def _to_class(self, cls, expr):

        if isinstance(expr, list):
            return ExprList(expr)
        elif isinstance(expr, tuple):
            return ExprTuple(expr)
        elif isinstance(expr, dict):
            return ExprDict(expr)
        return cls(expr)
    
# This will allow sym.sympify to magically extract the sympy expression
# but it will also bypass our __rmul__, __radd__, etc. methods that get called
# when sympy punts.  Thus pi * t becomes a Mul rather than tExpr.
#
#    def _sympy_(self):
#        # This is called from sym.sympify
#        return self.expr

    def __getattr__(self, attr):

        if False:
            print(self.__class__.__name__, attr)

        expr1 = self.expr            
        try:
            a = getattr(expr1, attr)
        except:
            raise

        # This gets called if there is no explicit attribute attr for
        # this instance.  We call the method of the wrapped sympy
        # class and rewrap the returned value if it is a sympy Expr
        # object.

        # FIXME.  This propagates the assumptions.  There is a
        # possibility that the operation may violate them.


        # If it is not callable, directly wrap it.
        if not callable(a):
            if not isinstance(a, sym.Expr):
                return a
            ret = a
            if hasattr(self, 'assumptions'):
                return self.__class__(ret, **self.assumptions)
            return self._to_class(self.__class__, ret)

        # If it is callable, create a function to pass arguments
        # through and wrap its return value.
        def wrap(*args, **kwargs):
            """This is quantity for a SymPy function.
            For help, see the SymPy documentation."""

            ret = a(*args, **kwargs)
            
            if not isinstance(ret, sym.Expr):
                # May have tuple, etc.   These could be wrapped but
                # it appears that this leads to more grief when working
                # with SymPy.
                return ret
            
            # Wrap the return value
            cls = self.__class__
            if hasattr(self, 'assumptions'):
                return cls(ret, **self.assumptions)
            return self._to_class(self.__class__, ret)            
        
        return wrap

    def debug(self):
        """Print the SymPy epxression and the assumptions for all symbols in
        the expression."""

        name = self.__class__.__name__
        s = '%s(' % name
        print(symdebug(self.expr, s, len(name) + 1))

    @property
    def expr_with_units(self):
        """Return SymPy expression with units."""

        if self.units == 1:
            return expr
        
        # Don't evaluate otherwise 1 A gets printed as A.
        return sym.Mul(self.expr, self.units, evaluate=False)

    @property
    def expr_with_canonical_units(self):
        """Return SymPy expression with canonical units."""

        return self.expr * self.canonical_units    
        
    @property
    def func(self):
        """Return the top-level function in the Sympy Expression.

        For example, this returns Mul for the expression `3 * s`.
        See also .args(), to return the args, in this case `(3, s)`"""

        return self.expr.func
    
    def __abs__(self):
        """Absolute value."""

        return self.__class__(self.abs, **self.assumptions)

    def __neg__(self):
        """Negation."""

        return self.__class__(-self.expr, **self.assumptions)

    def _incompatible(self, x, op, reason=''):

        raise ValueError("""
Cannot determine %s(%s) %s %s(%s)%s""" %
                         (self.__class__.__name__, self, op,
                          x.__class__.__name__, x, reason))        
                
    def _incompatible_domains(self, x, op):

        self._incompatible(x, op, ' since the domains are incompatible')

    def _incompatible_quantities(self, x, op):

        self._incompatible(x, op, """ since the units of the result are unsupported.
As a workaround use x.as_expr() %s y.as_expr()""" % op)

    def _add_compatible_domains(self, x):

        return self.domain == x.domain
    
    def _mul_compatible_domains(self, x):

        if self.domain == x.domain:
            return True
        
        if (self.is_constant_domain or x.is_constant_domain):
            return True

        return False

    def _div_compatible_domains(self, x):

        if self.domain == x.domain:
            return True
        
        if (self.is_constant_domain or x.is_constant_domain):
            return True

        return False

    def __compat_add__(self, x, op):

        assumptions = self.assumptions.copy()

        if not isinstance(x, Expr):
            x = expr(x)

        if state.check_units:
            sunits = self.canonical_units
            xunits = x.canonical_units
            
            if (sunits != xunits and self.expr != 0 and x.expr != 0 and not
                (state.loose_units and x.is_undefined)):
                self._incompatible(x, op, ' since the units %s are incompatible with %s' % (self.units, x.units))

        cls = self.__class__
        xcls = x.__class__
            
        if x.is_constant_domain and x.quantity == 'undefined':
            if state.loose_units or x.expr == 0:
                # Allow voltage(1) + 2 etc.
                return cls, self, x, assumptions
            if self.is_transfer:
                # Allow transfer(1) == 1
                return cls, self, x, assumptions                
        
        if self.is_constant_domain and self.quantity == 'undefined':
            return xcls, self, x, assumptions            

        if (isinstance(self, LaplaceDomainExpression) or
            isinstance(x, LaplaceDomainExpression)):
            assumptions = self.assumptions.add(x)
        
        if self.quantity == x.quantity:
            if self.is_constant_domain:
                return xcls, self, x, assumptions
            if x.is_constant_domain:            
                return cls, self, x, assumptions                
            if self.domain == x.domain:
                return cls, self, x, assumptions

        # For phasor comparisons...
        if self.is_phasor_domain and x.is_angular_fourier_domain:
            return cls, self, cls(x), assumptions
        if self.is_angular_fourier_domain and x.is_phasor_domain:
            return xcls, cls(self), x, assumptions        

        if not self._add_compatible_domains(x):
            self._incompatible_domains(x, op)

        # expr + voltage
        if self.quantity == 'undefined':
            if state.loose_units or x.is_transfer:
                return xcls, self, x, assumptions
            
        # voltage + expr
        if x.quantity == 'undefined':
            if state.loose_units or self.is_transfer:            
                return cls, self, x, assumptions        
        
        self._incompatible_quantities(x, op)        

    def __mul__(self, x):
        """Multiply."""

        from .super import Superposition

        if isinstance(x, Superposition):
            return x.__mul__(self)
        
        if not isinstance(x, Expr):
            x = expr(x)

        # Handle omega * t
        if (self.__class__ == AngularFourierDomainExpression and
            x.__class__ == TimeDomainExpression):
            return TimeDomainExpression(self.expr * x.expr)
        if (x.__class__ == AngularFourierDomainExpression and
            self.__class__ == TimeDomainExpression):
            return TimeDomainExpression(self.expr * x.expr)

        # Try to convert immitance to a constant so that can handle I(t) * Z
        if x.is_immitance:
            try:
                xunits = x.units
                x = x.as_constant()
                x.units = xunits
            except:
                pass

        if not self._mul_compatible_domains(x):
            self._incompatible_domains(x, '*')                    

        if self.is_transform_domain:
            assumptions = self.assumptions.convolve(x)
        else:
            assumptions = Assumptions()

        xquantity, yquantity = x.quantity, self.quantity
        # Maybe use undefined for voltage**2 etc.
        if xquantity == 'undefined':
            xquantity = 'constant'
        if yquantity == 'undefined':
            yquantity = 'constant'            
        
        key = (yquantity, xquantity)
        if key not in self._mul_mapping:
            key = (xquantity, yquantity)
            if key not in self._mul_mapping:
                # TODO: What about voltage**2. etc.
                self._incompatible_quantities(x, '*')        

        quantity = self._mul_mapping[key]
        if quantity == 'constant':
            quantity = 'undefined'

        if self.is_constant_domain:
            cls = x._class_by_quantity(quantity)
        else:
            cls = self._class_by_quantity(quantity)

        value = self.expr * x.expr
        result = cls(value, **assumptions)        
        result.units = self.units * x.units
        return result
    
    def __rmul__(self, x):
        """Reverse multiply."""

        return self.__mul__(x)

    def __truediv__(self, x):
        """True divide."""

        if not isinstance(x, Expr):
            x = expr(x)

        # Try to convert immitance to a constant so that can handle V(t) / Z
        if x.is_immitance:
            try:
                x = x.as_constant()
            except:
                pass
            
        if not self._div_compatible_domains(x):
            self._incompatible_domains(x, '/')                        

        assumptions = self.assumptions.convolve(x)

        xquantity, yquantity = x.quantity, self.quantity
        # Maybe use undefined for voltage**2 etc.
        if xquantity == 'undefined':
            xquantity = 'constant'
        if yquantity == 'undefined':
            yquantity = 'constant'            
        
        key = (yquantity, xquantity)
        if key not in self._div_mapping:
            # TODO: What about voltage**2. etc.
            self._incompatible_quantities(x, '/')        

        quantity = self._div_mapping[key]
        if quantity == 'constant':
            quantity = 'undefined'

        if self.is_constant_domain:
            cls = x._class_by_quantity(quantity)
        else:
            cls = self._class_by_quantity(quantity)

        value = self.expr / x.expr
        result = cls(value, **assumptions)        
        result.units = self.units / x.units

        return result
            
    def __rtruediv__(self, x):
        """Reverse true divide."""

        from .matrix import Matrix
        
        if isinstance(x, Matrix):
            return x / self.expr

        if not isinstance(x, Expr):
            x = expr(x)

        return x.__truediv__(self)
            
    def __add__(self, x):
        """Add."""

        from .matrix import Matrix

        if isinstance(x, Matrix):
            return x + self.expr

        cls, self, x, assumptions = self.__compat_add__(x, '+')
        return cls(self.expr + x.expr, **assumptions)

    def __radd__(self, x):
        """Reverse add."""

        if not isinstance(x, Expr):
            x = expr(x)
        return x.__add__(self)

    def __sub__(self, x):
        """Subtract."""

        from .super import Superposition

        if isinstance(x, Superposition):
            return -x + self
        
        cls, self, x, assumptions = self.__compat_add__(x, '-')
        return cls(self.expr - x.expr, **assumptions)

    def __rsub__(self, x):
        """Reverse subtract."""

        if not isinstance(x, Expr):
            x = expr(x)
        return x.__sub__(self)        

    def __pow__(self, x):
        """Power."""

        if x == 2:
            return self.__mul__(self)
        elif x == -1:
            return self.__rtruediv__(1)
        elif self.quantity != 'undefined':
            raise ValueError('Cannot compute %s(%s) ** %s' % (self.__class__.__name__, self, x))

        if not isinstance(x, Expr):
            x = expr(x)

        result = self.expr.__pow__(x.expr)
        if not self.is_constant_domain:
            return self.__class__(result)
        return x.__class__(result)        
    
    def __rpow__(self, x):
        """Reverse pow, x**self."""

        if not isinstance(x, Expr):
            x = expr(x)
        
        return x.__pow__(self)

    def __or__(self, x):
        """Parallel combination."""

        return self.parallel(x)

    def __eq__(self, x):
        """Test for mathematical equality as far as possible.
        This cannot be guaranteed since it depends on simplification.
        Note, SymPy comparison is for structural equality."""

        # Note, this is used by the in operator.

        if x is None:
            return False

        # Handle self == []
        if isinstance(x, list):
            return False

        # Disallow t == 't', etc.
        if isinstance(x, str):
            return False
        
        try:
            cls, self, x, assumptions = self.__compat_add__(x, '==')
        except ValueError:
            return False
            
        x = cls(x)

        # This fails if one of the operands has the is_real attribute
        # and the other doesn't...
        return sym.simplify(self.expr - x.expr) == 0

    def __ne__(self, x):
        """Test for mathematical inequality as far as possible.
        This cannot be guaranteed since it depends on simplification.
        Note, SymPy comparison is for structural equality."""

        if x is None:
            return True

        try:        
            cls, self, x, assumptions = self.__compat_add__(x, '!=')
        except ValueError:
            return True            
            
        x = cls(x)

        return sym.simplify(self.expr - x.expr) != 0        

    def __gt__(self, x):
        """Greater than."""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '>')
        x = cls(x)

        return self.expr > x.expr

    def __ge__(self, x):
        """Greater than or equal."""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '>=')
        x = cls(x)

        return self.expr >= x.expr

    def __lt__(self, x):
        """Less than."""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '<')
        x = cls(x)

        return self.expr < x.expr

    def __le__(self, x):
        """Less than or equal."""

        if x is None:
            return True

        cls, self, x, assumptions = self.__compat_add__(x, '<=')
        x = cls(x)

        return self.expr <= x.expr

    def parallel(self, x):
        """Parallel combination."""

        cls, self, x, assumptions = self.__compat_add__(x, '|')
        x = cls(x)

        return cls(self.expr * x.expr / (self.expr + x.expr), **assumptions)

    def copy(self):
        """Copy the expression."""
        return self.__class__(self.expr, **self.assumptions)

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
        if self.is_real:
            dst = self.__class__(0, **assumptions)
            dst.part = 'imaginary'
            return dst
        
        assumptions['real'] = True
        dst = self.__class__(symsimplify(sym.im(self.expr)), **assumptions)
        dst.part = 'imaginary'
        return dst

    @property
    def real_imag(self):
        """Rewrite as x + j * y."""

        return self.real + j * self.imag

    @property
    def _ratfun(self):

        try:
            return self.__ratfun
        except:
            pass

        if self.var is None:
            self.__ratfun = None
        else:
            # Note, this handles expressions that are products of
            # rational functions and arbitrary delays.  
            self.__ratfun = Ratfun(self.expr, self.var)
        return self.__ratfun

    @property
    def K(self):
        """Return gain."""

        return self.N.coeffs()[0] / self.D.coeffs()[0] 
    
    @property
    def N(self):
        """Return numerator of rational function.
        The denominator is chosen so that it is a polynomial."""

        return self.numerator

    @property
    def D(self):
        """Return denominator of rational function.
        The denominator is chosen so that it is a polynomial."""

        return self.denominator

    @property
    def numerator(self):
        """Return numerator of rational function.
        The denominator is chosen so that it is a polynomial."""

        N, D = self.as_N_D()
        return N

    @property
    def denominator(self):
        """Return denominator of rational function.
        The denominator is chosen so that it is a polynomial."""

        N, D = self.as_N_D()
        return D

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

    def divide_top_and_bottom(self, factor):
        """Divide numerator and denominator by common factor."""

        N = (self.N / factor).expand()
        D = (self.D / factor).expand()

        return N / D

    def factor_const(self):

        from .utils import factor_const

        c, r = factor_const(self, self.var)
        return ConstantDomainExpression(c), self.__class__(r, **self.assumptions)

    def term_const(self):

        from .utils import term_const

        c, r = term_const(self, self.var)
        return ConstantDomainExpression(c), self.__class__(r, **self.assumptions)    

    def multiply_top_and_bottom(self, factor):
        """Multiply numerator and denominator by common factor."""

        N = self.N.expr
        D = self.D.expr

        N = sym.Mul(N, factor, evaluate=False)
        D = sym.Mul(D, factor, evaluate=False)
        ID = sym.Pow(D, -1, evaluate=False)
        expr = sym.Mul(N, ID, evaluate=False)
        
        return self.__class__(expr)
    
    @property
    def magnitude(self):
        """Return magnitude."""

        if self.is_real:
            dst = expr(abs(self.expr))
            dst.part = 'magnitude'            
            return dst

        R = self.rationalize_denominator()
        N = R.N
        Dnew = R.D
        Nnew = sqrt((N.real**2 + N.imag**2).simplify())
        dst = Nnew / Dnew

        dst = dst.as_quantity(self.quantity)

        dst.part = 'magnitude'
        return dst

    @property
    def abs(self):
        """Return magnitude."""

        return self.magnitude

    @property
    def sign(self):
        """Return sign."""

        return self.__class__(sym.sign(self.expr), **self.assumptions)

    @property
    def dB(self):
        """Return magnitude in dB."""

        # Need to clip for a desired dynamic range?
        # Assume reference is 1.
        dst = 20 * log10(self.magnitude)
        dst.part = 'magnitude'
        dst.units = dB
        return dst

    @property
    def phase(self):
        """Return phase in radians."""

        if self.is_time_domain or self.is_discrete_time_domain:
            raise ValueError('Cannot determine phase of time-domain expression %s' % self)
        
        R = self.rationalize_denominator()
        N = R.N

        if N.imag == 0:
            if N.real >= 0:
                dst = expr(0)
            else:
                dst = expr(sym.pi)
        else:
            if N.real != 0:
                G = gcd(N.real, N.imag)
                N = N / G
            dst = atan2(N.imag, N.real)

        dst.part = 'phase'
        dst.units = uu.rad
        return dst

    @property
    def phase_degrees(self):
        """Return phase in degrees."""

        dst = self.phase * 180.0 / sym.pi
        dst.part = 'phase'
        dst.units = uu.deg
        return dst

    @property
    def angle(self):
        """Return phase angle (in radians)."""

        return self.phase

    @property
    def polar(self):
        """Return in polar format."""

        return self.abs * exp(j * self.phase)

    @property
    def cartesian(self):
        """Return in Cartesian format."""

        return self.real + j * self.imag
    
    @property
    def is_number(self):
        """Returns True if expression is a number."""                

        return self.expr.is_number

    @property
    def is_constant(self):
        """Returns True if expression does not have any free symbols  (compare with `is_unchanging`)."""        

        return self.expr.free_symbols == set()

    @property
    def is_unchanging(self):
        """Returns True if expression does not have a domain variable (compare with `is_constant`)."""

        return self.var not in self.expr.free_symbols
    
    def evaluate(self, arg=None):
        """Evaluate expression at arg.  arg may be a scalar, or a vector.
        The result is of type float or complex.

        There can be only one or fewer undefined variables in the expression.
        This is replaced by arg and then evaluated to obtain a result.
        """

        is_causal = self.is_causal
        
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

            def rect(arg):
                return 1.0 if abs(arg) <= 0.5 else 0.0

            def sinc(arg):
                return 1.0 if arg == 0 else np.sin(np.pi * arg) / (np.pi * arg)
            
            def dirac(arg):
                return np.inf if arg == 0.0 else 0.0

            def unitimpulse(arg):
                return 1.0 if arg == 0 else 0.0            

            def heaviside(arg):
                return 1.0 if arg >= 0.0 else 0.0

            def sqrt(arg):
                # Large numbers get converted to ints and int has no sqrt
                # attribute so convert to float.
                if isinstance(arg, int):
                    arg = float(arg)
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
            func1 = lambdify(var, expr,
                            ({'DiracDelta' : dirac,
                              'Heaviside' : heaviside,
                              'UnitImpulse' : unitimpulse,
                              'sinc' : sinc,
                              'rect' : rect,
                              'sqrt' : sqrt, 'exp' : exp},
                             "scipy", "numpy", "math", "sympy"))

            def func(arg):
                # Lambdify barfs on (-1)**n if for negative values of n.
                # even if have (-1)**n * Heaviside(n)
                # So this function heads Lambdify off at the pass,
                # if the function is causal.
                
                if is_causal and arg < 0:
                    return 0
                return func1(arg)

            try:
                result = func(arg0)
                response = complex(result)
            except NameError as e:
                raise RuntimeError('Cannot evaluate expression %s: %s' % (self, e))
            except AttributeError as e:
                if False and expr.is_Piecewise:
                    raise RuntimeError(
                        'Cannot evaluate expression %s,'
                        ' due to undetermined conditional result' % self)

                raise RuntimeError(
                    'Cannot evaluate expression %s,'
                    ' probably have a mysterious function: %s' % (self, e))

            except TypeError as e:
                raise RuntimeError('Cannot evaluate expression %s: %s' % (self, e))
            
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

        if not hasattr(self, 'var') or self.var is None:
            symbols = list(expr.free_symbols)
            if arg is None:
                if len(symbols) == 0:
                    return expr.evalf()
                raise ValueError('Undefined symbols %s in expression %s' % (tuple(symbols), self))                                    
            if len(symbols) == 0:
                print('Ignoring arg %s' % arg)
                return expr.evalf()
            elif len(symbols) == 1:            
                return evaluate_expr(expr, symbols[0], arg)
            else:
                raise ValueError('Undefined symbols %s in expression %s' % (tuple(symbols), self))                
                
            
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

        try:
            arg = arg.evalf()
        except:
            pass

        return evaluate_expr(expr, var, arg)

    def has(self, *patterns):
        """Test whether any subexpressions matches any of the patterns.  For example,
         V.has(exp(t)) 
         V.has(t)

        """

        tweak_patterns = [pattern.expr if isinstance(pattern, (Expr, Function)) else pattern for pattern in patterns]
        return self.expr.has(*tweak_patterns)

    def has_symbol(self, sym):
        """Test if have symbol contained.  For example,
        V.has_symbol('a')
        V.has_symbol(t)
        
        """                        
        return self.has(symbol(sym))
    
    def _subs1(self, old, new, **kwargs):

        # This will fail if a variable has different attributes,
        # such as positive or real.
        # Should check for bogus substitutions, such as t for s.

        if new is old:
            return self

        expr = new
        if isinstance(new, Expr):
            if old == self.var:
                domain = new.domain
            else:
                domain = self.domain
            expr = new.expr
        else:
            domain = self.domain
            expr = sympify(expr)

        old = symbol_map(old)

        if isinstance(expr, list):
            # Get lists from solve.  These stymy sympy's subs.
            if len(expr) == 1:
                expr = expr[0]
            else:
                print('Warning, substituting a list...')
        
        result = self.expr.subs(old, expr)

        # If get empty Piecewise, then result unknowable.  TODO: sympy
        # 1.2 requires Piecewise constructor to have at least one
        # pair.
        if False and result.is_Piecewise and result == sym.Piecewise():
            result = sym.nan

        return self.change(result, domain=domain, **self.assumptions)

    def transform(self, arg, **assumptions):
        """Transform into a different domain.

        If arg is f, s, t, omega, jomega perform domain transformation,
        otherwise perform substitution.

        Note (5 * s)(omega) will fail since 5 * s is assumed not to be
        causal and so Fourier transform is unknown.  However, Zs(5 *
        s)(omega) will work since Zs is assumed to be causal."""
        
        from .transform import transform
        return transform(self, arg, **assumptions)

    def __call__(self, arg, **assumptions):
        """Transform domain or substitute arg for variable. 
        
        Substitution is performed if arg is a tuple, list, numpy
        array, or constant.  If arg is a tuple or list return a list.
        If arg is an numpy array, return numpy array.

        Domain transformation is performed if arg is a domain variable
        or an expression of a domain variable.

        See also evaluate.

        """
        if isinstance(arg, (tuple, list)):
            return [self._subs1(self.var, arg1) for arg1 in arg]

        if isinstance(arg, np.ndarray):
            return np.array([self._subs1(self.var, arg1) for arg1 in arg])

        from .transform import call
        return call(self, arg, **assumptions)

    def _select(self, kind):

        from .transform import select
        return select(self, kind)        

    def limit(self, var, value, dir='+'):
        """Determine limit of expression(var) at var = value.
        If `dir == '+'` search from right else if `dir == '-'`
        search from left."""

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
        
        ret = sym.limit(self.expr, var, value, dir=dir)
        return self.__class__(ret, **self.assumptions)

    def simplify(self):
        """Simplify expression."""

        ret = symsimplify(self.expr)
        return self.__class__(ret, **self.assumptions)

    def simplify_units(self):
        """Simplify units into canonical form."""

        ret = self.__class__(self, **self.assumptions)
        ret.units = units.simplify_units(self.units)
        return ret

    def simplify_terms(self):
        """Simplify terms in expression individually."""

        result = 0
        for term in self.expr.expand().as_ordered_terms():
            result += symsimplify(term)
        return self.__class__(result, **self.assumptions)

    def simplify_factors(self):
        """Simplify factors in expression individually."""

        result = 0
        for factor in self.expr.factor().as_ordered_factors():
            result *= symsimplify(factor)
        return self.__class__(result, **self.assumptions)

    def simplify_sin_cos(self, as_cos=False, as_sin=False):
        """Simplify c * cos(theta) - s * sin(theta) as A * cos(theta - phi)."""

        result = simplify_sin_cos(self.expr, as_cos, as_sin)
        return self.__class__(result, **self.assumptions)        

    def replace(self, query, value, map=False, simultaneous=True, exact=None):

        try:
            query = query.expr
        except:
            pass

        try:
            value = value.expr
        except:
            pass        

        ret = self.expr.replace(query, value, map, simultaneous, exact)
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
        if hasattr(self, 'quantity_label'):
            label += self.quantity_label
            if self.part != '':
                label += ' ' + self.part
        else:
            if self.part != '':
                label += capitalize_name(self.part)
        return label

    @property
    def label_with_units(self):

        label = self.label
        if hasattr(self, 'units') and self.units != '' and self.units != 1:
            label += ' (%s)' % self.units
        return label    

    @property
    def domain_label_with_units(self):

        label = ''
        if hasattr(self, 'domain_label'):
            label += '%s' % self.domain_label
        if hasattr(self, 'domain_units'):
            if self.domain_units != '':
                label += ' (%s)' % self.domain_units
        return label

    def differentiate(self, arg=None):

        if arg is None:
            arg = self.var
        arg = self._tweak_arg(arg)
            
        return self.__class__(sym.diff(self.expr, arg), **self.assumptions)

    def diff(self, arg=None):

        return self.differentiate(arg)

    def _tweak_arg(self, arg):

        if isinstance(arg, Expr):
            return arg.expr

        if isinstance(arg, tuple):
            return tuple([self._tweak_arg(arg1) for arg1 in arg])

        if isinstance(arg, list):
            return [self._tweak_arg(arg1) for arg1 in arg]

        return arg

    def integrate(self, arg=None, **kwargs):

        if arg is None:
            arg = self.var

        arg = self._tweak_arg(arg)
        return self.__class__(sym.integrate(self.expr, arg, **kwargs),
                              **self.assumptions)

    def solve(self, *symbols, **flags):

        symbols = [symbol_map(symbol) for symbol in symbols]
        return expr(sym.solve(self.expr, *symbols, **flags))

    @property
    def symbols(self):
        """Return dictionary of symbols in the expression keyed by name."""
        symdict = {sym.name:sym for sym in self.free_symbols}

        # Look for V(s), etc.
        funcdict = {atom.func.__name__:atom for atom in self.atoms(AppliedUndef)}        

        symdict.update(funcdict)
        return symdict

    def _fmt_roots(self, roots, aslist=False):

        if not aslist:
            rootsdict = {}
            for root, n in roots.items():
                rootsdict[expr(root)] = n
            return expr(rootsdict)
            
        rootslist = []
        for root, n in roots.items():
            rootslist += [expr(root)] * n        
        return expr(rootslist)        
    
    def roots(self, aslist=False):
        """Return roots of expression as a dictionary
        Note this may not find them all."""

        if self._ratfun is None:
            roots = {}
        else:
            roots = self._ratfun.roots()
        return self._fmt_roots(roots, aslist)        
            
    def zeros(self, aslist=False):
        """Return zeroes of expression as a dictionary
        Note this may not find them all."""

        if self._ratfun is None:
            zeros = {}
        else:
            zeros = self._ratfun.zeros()
        return self._fmt_roots(zeros, aslist)        

    def poles(self, aslist=False, damping=None):
        """Return poles of expression as a dictionary
        Note this may not find them all."""

        if self._ratfun is None:
            return self._fmt_roots({}, aslist)            
        
        poles = self._ratfun.poles(damping=damping)

        polesdict = {}
        for pole in poles:
            key = pole.expr
            if key in polesdict:
                polesdict[key] += pole.n
            else:
                polesdict[key] = pole.n        

        return self._fmt_roots(polesdict, aslist)

    def canonical(self, factor_const=False):
        """Convert rational function to canonical form (aka polynomial form);
        this is like general form but with a unity highest power of
        denominator.  For example,

        (5 * s**2 + 5 * s + 5) / (s**2 + 4)

        If factor_const is True, factor constants from numerator, for example,

        5 * (s**2 + s + 1) / (s**2 + 4)

        This is also called gain-polynomial form.

        See also general, partfrac, standard, timeconst, and ZPK

        """
        if not self.expr.has(self.var):
            return self
        if self._ratfun is None:
            return self.copy()
        return self.__class__(self._ratfun.canonical(factor_const),
                              **self.assumptions)

    def general(self):
        """Convert rational function to general form.  For example,

        (5 * s**2 + 10 * s + 5) / (s**2 + 4)

        See also canonical, partfrac, standard, timeconst, and ZPK."""

        if self._ratfun is None:
            return self.copy()
        return self.__class__(self._ratfun.general(), **self.assumptions)

    def partfrac(self, combine_conjugates=False, damping=None):
        """Convert rational function into partial fraction form.   For example,

        5 + (5 - 15 * j / 4) / (s + 2 * j) + (5 + 15 * j / 4) / (s - 2 * j)

        If combine_conjugates is True then the pair of partial
        fractions for complex conjugate poles are combined.

        See also canonical, standard, general, timeconst, and ZPK."""

        try:
            if self._ratfun is None:
                return self.copy()        
            return self.__class__(self._ratfun.partfrac(combine_conjugates,
                                                        damping),
                                  **self.assumptions)
        except ValueError:
            return self.as_sum().partfrac(combine_conjugates, damping)

    def recippartfrac(self, combine_conjugates=False, damping=None):
        """Convert rational function into partial fraction form
        using reciprocal of variable.

        For example, if H = 5 * (s**2 + 1) / (s**2 + 5*s + 4)     
        then H.recippartfrac() gives 
        5/4 - 10/(3*(1 + 1/s)) + 85/(48*(1/4 + 1/s))

        If combine_conjugates is True then the pair of partial
        fractions for complex conjugate poles are combined.

        See also canonical, standard, general, partfrac, timeconst, and ZPK."""

        if self._ratfun is None:
            return self.copy()
        
        tmpsym = symsymbol('qtmp')

        expr = self.subs(1 / tmpsym)
        ratfun = Ratfun(expr.expr, tmpsym)

        nexpr = ratfun.partfrac(combine_conjugates, damping)
        nexpr = nexpr.subs(tmpsym, 1 / self.var)
        
        return self.__class__(nexpr, **self.assumptions)
    
    def standard(self):
        """Convert rational function into mixed fraction form.  For example,

        (5 * s - 5) / (s**2 + 4) + 5

        This is the sum of strictly proper rational function and a
        polynomial.

        See also canonical, general, partfrac, timeconst, and ZPK.

        """
        if self._ratfun is None:
            return self.copy()        
        return self.__class__(self._ratfun.standard(), **self.assumptions)

    def mixedfrac(self):
        """This is an alias for standard and may be deprecated."""
        
        return self.standard()

    def timeconst(self):
        """Convert rational function into time constant form.  For example,

        5 * (s**2 + 2 * s + 1) / (4 * (s**2 / 4 + 1))

        See also timeconst_terms, canonical, general, standard,
        partfrac and ZPK."""

        if self._ratfun is None:
            return self.copy()        
        return self.__class__(self._ratfun.timeconst(), **self.assumptions)

    def timeconst_terms(self):
        """Convert each term of expression into time constant form."""

        result = 0
        for term in self.expr.as_ordered_terms():
            result += self.__class__(term).timeconst()
        return self.__class__(result, **self.assumptions)            

    def ZPK(self):
        """Convert to zero-pole-gain (ZPK) form (factored form).  For example,

        5 * (s + 1)**2 / ((s - 2 * j) * (s + 2 * j))

        Note, both the numerator and denominator are expressed as
        products of monic factors, i.e., (s + 1 / 3) rather than (3 * s + 1).

        See also canonical, general, standard, partfrac, and timeconst.

        """

        if self._ratfun is None:
            return self.copy()        
        return self.__class__(self._ratfun.ZPK(), **self.assumptions)

    def factored(self):
        """Convert to factored form.  For example,

        5 * (s + 1)**2 / ((s - 2 * j) * (s + 2 * j))

        This is an alias for ZPK.  See also canonical, general,
        standard, partfrac, and timeconst.

        """
        
        if self._ratfun is None:
            return self.copy()
        return self.__class__(self._ratfun.ZPK(), **self.assumptions)
    
    def expandcanonical(self):
        """Expand in terms for different powers with each term
        expressed in canonical form.  For example,

        s / (s**2 + 4) + 5 / (s**2 + 4)

        See also canonical, general, partfrac, timeconst, and ZPK."""

        if self._ratfun is None:
            return self.copy()        
        return self.__class__(self._ratfun.expandcanonical(), **self.assumptions)

    def coeffs(self, norm=False):
        """Return list of coeffs assuming the expr is a polynomial in s.  The
        highest powers come first.  This will fail for a rational function.
        Instead use expr.N.coeffs or expr.D.coeffs for numerator
        or denominator respectively.
        
        If norm is True, normalise coefficients to highest power is 1."""

        if self._ratfun is None:
            return expr([self])
        
        try:
            z = sym.Poly(self.expr, self.var)
        except:
            raise ValueError('Use .N or .D attribute to specify numerator or denominator of rational function')

        c = z.all_coeffs()
        if norm:
            return expr([sym.simplify(c1 / c[0]) for c1 in c])
            
        return expr(c)

    def normcoeffs(self):
        """Return list of coeffs (normalised so the highest power is 1)
        assuming the expr is a polynomial in s.  The highest powers
        come first.  This will fail for a rational function.  Instead
        use expr.N.normcoeffs or expr.D.normcoeffs for numerator or
        denominator respectively."""

        return self.coeffs(norm=True)

    @property
    def degree(self):
        """Return the degree (order) of the rational function.

        This the maximum of the numerator and denominator degrees.
        Note zero has a degree of -inf."""

        if self._ratfun is None:
            return 1
        
        return self._ratfun.degree

    @property
    def Ndegree(self):
        """Return the degree (order) of the numerator of a rational function.
        This will throw an exception if the expression is not a
        rational function.

        Note zero has a degree of -inf.

        """

        if self._ratfun is None:
            return 1
        
        return self._ratfun.Ndegree

    @property
    def Ddegree(self):
        """Return the degree (order) of the denominator of a rational function.
        This will throw an exception if the expression is not a
        rational function.

        Note zero has a degree of -inf."""

        if self._ratfun is None:
            return 1        
        
        return self._ratfun.Ddegree

    def prune_HOT(self, degree):
        """Prune higher order terms if expression is a polynomial
        so that resultant approximate expression has the desired degree."""

        coeffs = self.coeffs
        if len(coeffs) < degree:
            return self

        coeffs = coeffs[::-1]

        expr = sym.S.Zero
        var = self.var
        for m in range(degree + 1):
            term = coeffs[m].expr * var ** m
            expr += term

        return self.__class__(expr, **self.assumptions)            

    def ratfloat(self):
        """This converts rational numbers in an expression to floats.
        See also floatrat.

        For example, t / 5 -> 0.2 * t
        """

        expr = self.expr
        expr = expr.replace(lambda expr: expr.is_Rational,
                            lambda expr: sym.Float(expr))

        return self.__class__(expr, **self.assumptions)

    def floatrat(self):
        """This converts floating point numbers to rational numbers in an
        expression.  See also ratfloat.

        For example, 0.2 * t - > t / 5

        """

        expr = self.expr
        # FIXME...
        expr = expr.replace(lambda expr: expr.is_Float,
                            lambda expr: sym.Rational(expr))

        return self.__class__(expr, **self.assumptions)            
    
    def approximate_fractional_power(self, order=2):
        """This is an experimental method to approximate
        s**a, where a is fractional, with a rational function using
        a Pade approximant."""

        v = self.var
        
        def query(expr):

            if not expr.is_Pow:
                return False
            if expr.args[0] != v:
                return False
            if expr.args[1].is_Number and not expr.args[1].is_Integer:
                return True
            if expr.args[1].is_Symbol and not expr.args[1].is_Integer:
                return True
            return False

        def value1(expr):

            a = expr.args[1]

            n = v * (a + 1) + (1 - a)
            d = v * (a - 1) + (1 + a)
            return n / d

        def value2(expr):

            a = expr.args[1]

            n = v**2 * (a**2 + 3 * a + 2) + v * (8 - a**2) + (a**2 - 3 * a + 2)
            d = v**2 * (a**2 - 3 * a + 2) + v * (8 - a**2) + (a**2 + 3 * a + 2)
            return n / d        

        if order == 1:
            value = value1
        elif order == 2:
            value = value2
        else:
            raise ValueError('Can only handle order 1 and 2 at the moment')
        
        expr = self.expr
        expr = expr.replace(query, value)

        return self.__class__(expr, **self.assumptions)

    def as_value_unit(self):
        """Return tuple of value and unit.  For example,

        >>> v = voltage(5)
        >>> v.as_value_unit
        (5, volts)
        """
        
        from .units import units

        return units.as_value_unit(self.expr)
            
    def as_N_D(self, monic_denominator=False):
        """Responses due to a sum of delayed transient responses
        cannot be factored into ZPK form with a constant delay.
        For example, sometimes SymPy gives:

            ⎛    s⋅τ     ⎞  -s⋅τ
            ⎝V₁⋅ℯ    - V₂⎠⋅ℯ    
        I = ────────────────────
               s⋅(L⋅s + R)     

        This method tries to extract the numerator and denominator
        where the denominator is a polynomial.

        N, D = I.as_N_D()

                     -s⋅τ
        N = V₁ - V₂⋅ℯ    
        D =  s⋅(L⋅s + R)"""

        N, D = as_N_D(self.expr, self.var, monic_denominator)

        # Strip quantity and assumptions
        cls = self._class_by_quantity('undefined')
        return cls(N), cls(D)

    def as_sum(self):
        """Responses due to a sum of delayed transient responses
        cannot be factored into ZPK form with a constant delay.
        For example, sometimes SymPy gives:

            ⎛    s⋅τ     ⎞  -s⋅τ
            ⎝V₁⋅ℯ    - V₂⎠⋅ℯ    
        I = ────────────────────
               s⋅(L⋅s + R)     

        While this cannot be factored into ZPK form, it can be
        expressed as a sum of ZPK forms or as a partial fraction
        expansion.  However, SymPy does not play ball if trying to
        express as a sum of terms:

        I.as_ordered_terms()  
                                                 
        ⎡⎛    s⋅τ     ⎞  -s⋅τ⎤
        ⎢⎝V₁⋅ℯ    - V₂⎠⋅ℯ    ⎥
        ⎢────────────────────⎥
        ⎣    s⋅(L⋅s + R)     ⎦

        Instead, it appears necessary to split into N / D where
        D is a polynomial.  Then N can be split.
        """

        result = as_sum(self.expr, self.var)
        return self.__class__(result, **self.assumptions)

    def as_monic_terms(self):
        """Rewrite terms so that each denominator is monic.

        This does not expand the expression first; use `.expand()`."""

        result = 0
        for term in self.expr.as_ordered_terms():
            N, D = as_N_D(term, self.var, monic_denominator=True)
            result += N / D
        return self.__class__(result, **self.assumptions)

    def as_nonmonic_terms(self):
        """Rewrite terms so that each denominator is not monic.

        This does not expand the expression first; use `.expand()`."""

        result = 0
        for term in self.expr.as_ordered_terms():
            N, D = as_N_D(term, self.var, monic_denominator=False)
            result += N / D
        return self.__class__(result, **self.assumptions)                

    def continued_fraction_coeffs(self):

        coeffs = []
        var = self.var        
        
        def foo(Npoly, Dpoly):

            # This seems rather complicated to extract the leading terms.
            NLM, NLC = Npoly.LT()
            DLM, DLC = Dpoly.LT()
            NLT = sym.Poly(NLM.as_expr() * NLC, var)
            DLT = sym.Poly(DLM.as_expr() * DLC, var)

            Q = NLT / DLT
            coeffs.append(Q)

            Npoly2 = sym.Poly(Npoly.as_expr() - Q * Dpoly.as_expr(), var)            
            if Npoly2 != 0:
                foo(Dpoly, Npoly2)

        N, D = self.expr.as_numer_denom()
        Npoly = sym.Poly(N, var)
        Dpoly = sym.Poly(D, var)

        if Dpoly.degree() > Npoly.degree():
            coeffs.append(0)
            Npoly, Dpoly = Dpoly, Npoly
        
        foo(Npoly, Dpoly)

        return expr(coeffs)
    
    def as_continued_fraction(self):
        """Convert expression into acontinued fraction."""

        def foo(coeffs):

            if len(coeffs) == 1:
                return coeffs[0]
            return coeffs[0] + 1 / foo(coeffs[1:])

        coeffs = self.continued_fraction_coeffs()
        result = foo(coeffs)
        return self.__class__(result, **self.assumptions)

    def continued_fraction_inverse_coeffs(self):

        coeffs = []
        var = self.var
        
        def foo(Npoly, Dpoly):

            # This seems rather complicated to extract the last non-zero terms.
            NEM, NEC = Npoly.ET()
            DEM, DEC = Dpoly.ET()
            NET = NEM.as_expr() * NEC
            DET = DEM.as_expr() * DEC
            
            if sym.Poly(NET, var).degree() > sym.Poly(DET, var).degree():
                coeffs.append(0)
                foo(Dpoly, Npoly)
                return

            Q = NET / DET
            coeffs.append(Q)
            
            Npoly2 = sym.Poly(Npoly.as_expr() - Q * Dpoly.as_expr(), var)
            if Npoly2 != 0:
                foo(Dpoly, Npoly2)

        N, D = self.expr.as_numer_denom()
        Npoly = sym.Poly(N, var)
        Dpoly = sym.Poly(D, var)

        foo(Npoly, Dpoly)
        return expr(coeffs)

    def as_continued_fraction_inverse(self):

        def foo(coeffs):

            if len(coeffs) == 1:
                return coeffs[0]
            return coeffs[0] + 1 / foo(coeffs[1:])

        coeffs = self.continued_fraction_inverse_coeffs()
        result = foo(coeffs)
        return self.__class__(result, **self.assumptions)

    def change(self, expr, domain=None, units_scale=None, **assumptions):
        """Change expression class."""

        if domain is None:
            domain = self.domain

        quantity = self.quantity

        cls = self._class_by_quantity(quantity, domain)        
        ret = cls(expr, **assumptions)

        if units_scale is not None:
            ret.units = self.units * units_scale
        return ret

    def force_time(self):

        if state.force_time:
            return self.time()
        return self

def exprcontainer(arg, **assumptions):

    if isinstance(arg, (ExprList, ExprTuple, ExprDict)):
        return arg
    elif isinstance(arg, list):
        return ExprList(arg)
    elif isinstance(arg, tuple):
        return ExprTuple(arg)
    elif isinstance(arg, dict):
        return ExprDict(arg)    
    elif isinstance(arg, np.ndarray):
        from .vector import Vector
        if arg.ndim > 1:
            raise ValueError('Multidimensional arrays unsupported; convert to Matrix')
        return Vector(arg)
    
    raise ValueError('Unsupported exprcontainer %s' % arg.__class__.name)

    
def expr(arg, **assumptions):
    """Create Lcapy expression from arg.

    If `arg` is an `Expr` it is returned, unless `assumptions` is specified.

    If `arg` is a string:
       If a t symbol is found in the string a TimeDomainExpression object is created.
       If a s symbol is found in the string a LaplaceDomainExpression object is created.
       If a f symbol is found in the string an FourierDomainExpression object is created.
       If an omega symbol is found in the string an AngularFourierDomainExpression object is created.

    For example, v = expr('3 * exp(-t / tau) * u(t)')

    V = expr('5 * s', causal=True)
    """

    from .sym import tsym, fsym, ssym, omegasym

    if arg is None:
        return arg
    
    if isinstance(arg, Expr) and assumptions == {}:
        return arg
    
    if not isinstance(arg, str) and hasattr(arg, '__iter__'):
        return exprcontainer(arg)
    
    expr = sympify(arg, **assumptions)

    symbols = expr.free_symbols
    
    if tsym in symbols:
        return texpr(expr, **assumptions)
    elif ssym in symbols:
        return sexpr(expr, **assumptions)
    elif fsym in symbols:
        return fexpr(expr, **assumptions)
    elif omegasym in symbols:
        return omegaexpr(expr, **assumptions)
    else:
        return cexpr(expr, **assumptions)


def symbol(name, **assumptions):
    """Create an Lcapy symbol.

    By default, symbols are assumed to be positive unless real is
    defined or positive is defined as False."""
    return expr(symsymbol(name, **assumptions))


def symbols(names, **assumptions):
    """Create Lcapy symbols from whitespace or comma delimiter string of
    symbol names.  See also symbol."""

    from .parser import split

    namelist = split(names, ", ")
    symbols = []
    for name in namelist:
        symbols.append(symbol(name, **assumptions))
    return symbols


def radians(arg, **assumptions):
    """Convert degrees to radians and set units to radians.  See also rad()
    that sets units as radians."""

    expr1 = expr(arg / 180 * sym.pi, **assumptions)
    expr1.units = uu.rad
    return expr1


def degrees(arg, **assumptions):
    """Convert radians to degrees and set units to degrees.  See also deg()
    that sets units as degrees."""

    expr1 = expr(arg * 180 / sym.pi, **assumptions)
    expr1.units = uu.deg
    return expr1


def rad(arg, **assumptions):
    """Set units to radians.  See also radians() that converts degrees to radians."""

    expr1 = expr(arg, **assumptions)
    expr1.units = uu.rad
    return expr1


def deg(arg, **assumptions):
    """Set units to degrees.  See also degrees() that converts radians to degrees."""

    expr1 = expr(arg, **assumptions)
    expr1.units = uu.deg
    return expr1


from .cexpr import cexpr, ConstantDomainExpression
from .fexpr import f, fexpr, FourierDomainExpression
from .texpr import t, texpr, TimeDomainExpression
from .sexpr import s, sexpr, LaplaceDomainExpression
from .omegaexpr import omega, omegaexpr, AngularFourierDomainExpression
from .expressionclasses import expressionclasses

# Horrible hack to work with IPython around Sympy's back for LaTeX
# formatting.  The problem is that Sympy does not check for the
# _repr_latex method and instead relies on a predefined list of known
# types.  See _can_print_latex method in sympy/interactive/printing.py

import sys
try:
    from .printing import latex
    formatter = sys.displayhook.shell.display_formatter.formatters['text/latex']
    
    for cls in (ExprList, ExprTuple, ExprDict):
        formatter.type_printers[cls] = Expr._repr_latex_
except:
    pass
