# -*- coding: utf-8 -*-

"""This module provides the Expr class.  This attempts to create a
consistent interface to SymPy's expressions.

Copyright 2014--2022 Michael Hayes, UCECE

"""


# TODO, propagate assumptions for arithmetic.........  This may be
# tricky.  At the moment only a limited propagation of assumptions are
# performed.

from __future__ import division
import sys
from .assumptions import Assumptions
from .cache import cached_property
from .domains import UndefinedDomain
from .quantity import UndefinedQuantity
from .ratfun import Ratfun
from .sym import sympify, symsimplify, j, omegasym, symdebug, AppliedUndef
from .sym import capitalize_name, tsym, miscsymbol, usersymbol, symbol_map, tausym, nusym, oo
from .sym import fsym, ssym, Fsym, Omegasym, symbol_delete, pi
from .sym import nsym, ksym, zsym
from .state import state, validate
from .printing import pprint, pretty, print_str, latex
from .functions import sqrt, log10, atan2, gcd, exp, Function, Eq
from .units import units, u as uu, dB
from .utils import (as_N_D, as_sum, remove_images, pair_conjugates,
                    split_dirac_delta, expand_functions)
import sympy as sym
from sympy.utilities.lambdify import lambdify
from sympy import UnevaluatedExpr
from .sym import simplify
from .simplify import simplify_sin_cos, simplify_heaviside, simplify_dirac_delta
from .simplify import simplify_rect, simplify_unit_impulse, simplify_conjugates
from .simplify import expand_hyperbolic_trig
from .approximate import (approximate_fractional_power, approximate_exp,
                          approximate_hyperbolic_trig, approximate_dominant,
                          approximate_order)
from .config import heaviside_zero, unitstep_zero
from collections import OrderedDict
from warnings import warn


__all__ = ('expr', 'symbol', 'symbols', 'deg', 'rad',
           'equation', 'UnevaluatedExpr')


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

    def _repr_png_(self):

        # This adds a large amount of data to the notebook without
        # adding anything useful.  It also triggers a bug on Windows
        # when converting a latex expression to png.
        return None

    def _repr_svg_(self):

        return None

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

    def latex_with_units(self, eng_format=False, show_units=True,
                         evalf=True, num_digits=3, **kwargs):
        """Make latex string with optional units."""

        from .engformatter import EngFormatter

        expr = self

        if evalf:
            expr = expr.evalf(num_digits)

        if show_units:
            units = str(expr.units)
            if units == '1':
                units = ''
        else:
            units = ''

        value = expr.sympy

        if evalf and value.is_number and eng_format:
            return EngFormatter(num_digits=num_digits).latex(value, units)

        s = latex(value, **kwargs)
        if show_units and units != '':
            s += '\,' + units
        return s

    def latex_math(self, **kwargs):
        """Make latex math-mode string."""
        return '$' + self.latex(**kwargs) + '$'

    def latexans(self, name, **kwargs):
        """Print latex string with LHS name."""
        return latex(sym.Eq(sympify(name), self._pexpr), **kwargs)

    def srepr(self):
        return sym.repr(self)


class ExprContainer(object):

    @property
    def sympy(self):
        """Return SymPy expression."""

        return self.expr

    def evaluate(self):
        """Evaluate each element to convert to floating point.
        This may change..."""

        return self.__class__([v.evalf() for v in self])

    def evalf(self, n=15):
        """Evaluate each element to convert to floating point values.
        `n` is the number of decimal places."""

        return self.__class__([v.evalf(n) for v in self])

    def simplify(self, **kwargs):
        """Simplify each element."""

        return self.__class__([simplify(v, **kwargs) for v in self])

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

        import pdb
        pdb.set_trace()
        return self


class ExprDict(ExprPrint, ExprContainer, ExprMisc, OrderedDict):
    """Facade class for dictionary created by sympy."""

    def __getitem__(self, key):

        # This is used for nodalanalysis to store results
        # indexed by node name.
        key2 = str(key)
        try:
            return super(ExprDict, self).__getitem__(key2)
        except:
            return super(ExprDict, self).__getitem__(key)

    def __call__(self, *args, **kwargs):
        """Perform substitution/transformation on each element."""

        new = {}
        for key, val in self.items():
            new[key] = expr(val)(*args, **kwargs)
        return self.__class__(new)

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

    def simplify(self, **kwargs):
        """Simplify each element but not the keys."""

        new = self.__class__()
        for k, v in self.items():
            new[k] = simplify(v, **kwargs)
        return new

    def has(self, *patterns):
        """Test whether any subexpressions matches any of the patterns.  For example,
         V.has(exp(t))
         V.has(t)

        """

        for value in self.values():
            if value.has(*patterns):
                return True
        return False

    def remove_undefs(self, return_mappings=False):
        """Remove applied undefined functions (undefs) with symbols.
        For example, `V(s)` is replaced with `V`."""

        mappings = {}
        new = {}
        for key, value in self.items():
            new[key], mappings1 = expr(
                value).remove_undefs(return_mappings=True)
            mappings.update(mappings1)

        new = self.__class__(new)
        if return_mappings:
            return new, mappings
        else:
            return new

    def solve(self, *symbols, **flags):
        """Solve system of equations, and return as ExprDict.
        See sympy.solve for usage."""

        if self.has(AppliedUndef):
            new, defs = self.remove_undefs(return_mappings=True)
            return new.solve(*symbols, **flags).subs(defs)

        symbols = delcapify(ExprTuple(symbols).remove_undefs())
        system = list(delcapify(self).values())
        solutions = sym.solve(system, *symbols, **flags)
        return expr(solutions)

    def evalf(self, n=15):
        """Evaluate each element to convert to floating point values.
        `n` is the number of decimal places."""

        new = self.__class__()
        for k, v in self.items():
            try:
                k = k.evalf(n)
            except:
                pass
            try:
                v = v.evalf(n)
            except:
                pass
            new[k] = v
        return new

    def subs(self, *args, **kwargs):
        """Substitute symbols for each expression in dict.

        `args` is either:
        - two arguments, e.g. foo.subs(old, new)
        - one iterable argument, e.g. foo.subs(iterable). The iterable may be
           o an iterable container with (old, new) pairs. In this case the
             replacements are processed in the order given with successive
             patterns possibly affecting replacements already made.
           o a dict or set whose key/value items correspond to old/new pairs.

        The domain of the result can be coerced using `domain` argument."""

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

    @property
    def expr(self):

        # The behaviour of this may change.  Perhaps the keys
        # should stay the same?
        new = {}
        for k, v in self.items():
            if isinstance(k, Expr):
                k = k.expr
            if isinstance(v, Expr):
                v = v.expr
            new[k] = v
        return new


class ExprList(ExprPrint, list, ExprContainer, ExprMisc):
    """Facade class for list created by sympy."""

    # Have ExprPrint first so that its _repr__pretty_ is called
    # in preference to list's one.  Alternatively, add explicit
    # _repr_pretty_ method here.

    def __init__(self, iterable=None, evalf=False, **assumptions):

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
                item = expr(item, **assumptions)
            eiterable.append(item)

        super(ExprList, self).__init__(eiterable)

    def __call__(self, *args, **kwargs):
        """Perform substitution/transformation on each element."""

        ret = [elt(*args, **kwargs) for elt in self]
        return self.__class__(ret)

    def subs(self, *args, **kwargs):
        """Substitute symbols for each expression in list.

        `args` is either:
        - two arguments, e.g. foo.subs(old, new)
        - one iterable argument, e.g. foo.subs(iterable). The iterable may be
           o an iterable container with (old, new) pairs. In this case the
             replacements are processed in the order given with successive
             patterns possibly affecting replacements already made.
           o a dict or set whose key/value items correspond to old/new pairs.

        The domain of the result can be coerced using `domain` argument."""

        return expr([e.subs(*args, **kwargs) for e in self])

    def has(self, *patterns):
        """Test whether any subexpressions matches any of the patterns.  For example,
         V.has(exp(t))
         V.has(t)

        """

        for value in self:
            if value.has(*patterns):
                return True
        return False

    def remove_undefs(self, return_mappings=False):
        """Remove applied undefined functions (undefs) with symbols.
        For example, `V(s)` is replaced with `V`."""

        mappings = {}
        new = []
        for value in self:
            new1, mappings1 = value.remove_undefs(return_mappings=True)
            new.append(new1)
            mappings.update(mappings1)

        new = self.__class__(new)
        if return_mappings:
            return new, mappings
        else:
            return new

    def solve(self, *symbols, **flags):
        """Solve system of equations and return as ExprDict.
        See sympy.solve for usage."""

        if self.has(AppliedUndef):
            new, defs = self.remove_undefs(return_mappings=True)
            return new.solve(*symbols, **flags).subs(defs)

        symbols = delcapify(ExprTuple(symbols).remove_undefs())
        system = delcapify(self)
        solutions = sym.solve(system, *symbols, **flags)
        return expr(solutions)

    @property
    def expr(self):
        return [e.expr for e in self]

    @property
    def fval(self):
        """Evaluate expression and return as a list of python float values."""

        return [float(a.fval) for a in self]

    @property
    def cval(self):
        """Evaluate expression and return as a list of python complex values."""

        return [complex(a.cval) for a in self]


class ExprTuple(ExprPrint, tuple, ExprContainer, ExprMisc):
    """Facade class for tuple created by sympy."""

    # Tuples are immutable, need to use __new__
    def __new__(cls, iterable, **assumptions):

        eiterable = [expr(e, **assumptions) for e in iterable]
        return super(ExprTuple, cls).__new__(cls, eiterable)

    def __call__(self, *args, **kwargs):
        """Perform substitution/transformation on each element."""

        ret = [elt(*args, **kwargs) for elt in self]
        return self.__class__(ret)

    def subs(self, *args, **kwargs):
        """Substitute symbols for each expression in tuple.

        `args` is either:
        - two arguments, e.g. foo.subs(old, new)
        - one iterable argument, e.g. foo.subs(iterable). The iterable may be
           o an iterable container with (old, new) pairs. In this case the
             replacements are processed in the order given with successive
             patterns possibly affecting replacements already made.
           o a dict or set whose key/value items correspond to old/new pairs.

        The domain of the result can be coerced using `domain` argument."""

        return expr(tuple([e.subs(*args, **kwargs) for e in self]))

    def has(self, *patterns):
        """Test whether any subexpressions matches any of the patterns.  For example,
         V.has(exp(t))
         V.has(t)

        """

        for value in self:
            if value.has(*patterns):
                return True
        return False

    def remove_undefs(self, return_mappings=False):
        """Remove applied undefined functions (undefs) with symbols.
        For example, `V(s)` is replaced with `V`."""

        mappings = {}
        new = ()
        for value in self:
            new1, mappings1 = value.remove_undefs(return_mappings=True)
            new = new + (new1, )
            mappings.update(mappings1)

        new = self.__class__(new)
        if return_mappings:
            return new, mappings
        else:
            return new

    def solve(self, *symbols, **flags):
        """Solve system of equations, and return as ExprDict.
        See sympy.solve for usage."""

        if self.has(AppliedUndef):
            new, defs = self.remove_undefs(return_mappings=True)
            return new.solve(*symbols, **flags).subs(defs)

        symbols = delcapify(ExprTuple(symbols).remove_undefs())
        system = delcapify(self)
        solutions = sym.solve(system, *symbols, **flags)
        return expr(solutions)

    @property
    def expr(self):
        return tuple([e.expr for e in self])

    @property
    def fval(self):
        """Evaluate expression and return as a tuple of python float values."""

        return tuple([float(a.fval) for a in self])

    @property
    def cval(self):
        """Evaluate expression and return as a tuple of python complex values."""

        return tuple([complex(a.cval) for a in self])


class ExprDomain(object):

    is_sequence = False

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
        elif domain == 'frequency response':
            return self.as_frequency_response()
        elif domain == 'angular frequency response':
            return self.as_angular_frequency_response()
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
        raise ValueError(
            'Can only convert voltage or current to superposition')

    def change(self, arg, domain=None, units_scale=None, **assumptions):
        """Change expression class."""

        if domain is None:
            domain = self.domain

        if domain == 'constant':
            # Allow changing of constants, e.g., V1 to 5 * t
            domain = expr(arg).domain

        quantity = self.quantity

        cls = self._class_by_quantity(quantity, domain)
        ret = cls(arg, **assumptions)

        if units_scale is not None:
            ret.units = self.units * units_scale
        return ret


class Expr(UndefinedDomain, UndefinedQuantity, ExprPrint, ExprMisc, ExprDomain):
    """Facade class for sympy classes derived from sympy.Expr."""

    var = None

    # This provides a minor speed improvement for attribute lookup
    # but prevents dynamic adding of new attributes.
    __slots__ = ()

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

    def __init__(self, arg, rational=True, **assumptions):
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
        # if arg == 0:
        #    assumptions.set('causal', True)
        if self.is_always_causal:
            assumptions.set('causal', True)

        self.assumptions = assumptions
        # Remove Lcapy assumptions from SymPy expr.
        self.expr = sympify(arg, rational=rational, **
                            self.assumptions.sympy_assumptions())
        try:
            self._units = self._default_units
        except:
            self._units = sym.S.One

    def as_time(self):
        return self.time()

    def as_laplace(self):
        return self._cached_laplace()

    def as_phasor(self):
        return self.phasor()

    def as_fourier(self):
        return self.fourier()

    def as_angular_fourier(self):
        return self.angular_fourier()

    def as_frequency_response(self):
        return self.frequency_response()

    def as_angular_frequency_response(self):
        return self.angular_frequency_response()

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
        s = latex(self._pexpr, mode='plain')
        return "$$%s$$" % s

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
        """Return the canonical units of the expression.  This is a simplified
        form, so volt * ampere becomes watt.

        """

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
    def is_complex_signal(self):
        """Return True if time-domain signal is complex."""

        if 'complex_signal' not in self.assumptions:
            return False
        return self.assumptions['complex_signal'] == True

    @property
    def is_complex(self):
        from .sym import ssym
        from .sym import zsym

        if self.has(ssym) or self.has(zsym):
            return True

        # Sometimes there is a lingering Re or Im operator
        # even though we know the result is real.
        if self.part != '':
            return False

        return self.has(j)

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
        """Evaluate expression and return as a python float value.
        Small imaginary components are discarded.

        Use `cval` to return a complex value.
        """

        val = self.val.expr
        try:
            return float(val)
        except TypeError:
            pass
        real, imag = val.as_real_imag()
        real = float(real)
        imag = float(imag)
        if abs(imag / real) > 1e-15:
            warn('Discarding imaginary part')
        return real

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

        new = self.copy()
        # Don't create Expr since SymPy sympify will create Integers
        # rather than Floats if the truncated Float looks like an integer.
        new.expr = self.ratfloat().expr.evalf(n, *args, **kwargs)
        return new

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
# when sympy punts.  Thus pi * t becomes a Mul rather than TimeDomainExpression.
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
            # Hack for ubuntu-20.04, python 3.7 and 3.8
            if attr == 'abbrev':
                return ''
            raise AttributeError("'%s' object has no attribute '%s'" % (
                self.__class__.__name__, attr))

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

            # Extract SymPy expressions from Lcapy expressions
            newargs = []
            for arg in args:
                try:
                    newargs.append(arg.expr)
                except AttributeError:
                    newargs.append(arg)

            # Extract SymPy expressions from Lcapy expressions
            newkwargs = {}
            for key, arg in kwargs.items():
                try:
                    newkwargs[key] = kwargs[key].expr
                except AttributeError:
                    newkwargs[key] = kwargs[key]

            ret = a(*newargs, **newkwargs)

            if not isinstance(ret, sym.Expr):
                # Hack for jupyter notebook printer returning
                # png as a bytes object
                if isinstance(ret, bytes):
                    return ret
                # Wrap list, tuple, etc.
                return expr(ret)

            # Wrap the return value
            cls = self.__class__
            if hasattr(self, 'assumptions'):
                return cls(ret, **self.assumptions)
            return self._to_class(self.__class__, ret)

        return wrap

    def debug(self):
        """Print the SymPy expression and the assumptions for all symbols in
        the expression.  See also srepr."""

        name = self.__class__.__name__
        s = '%s(' % name
        print(symdebug(self.expr, s, len(name) + 1))

    def srepr(self):
        """Print the SymPy abstract syntax tree for the expression."""

        sym.srepr(self.sympy)

    @property
    def sympy(self):
        """Return SymPy expression."""

        return self.expr

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
        See also args to return the args, in this case `(3, s)`"""

        return self.expr.func

    def __abs__(self):
        """Absolute value."""

        return self.__class__(self.abs, **self.assumptions)

    def __neg__(self):
        """Negation."""

        return self.__class__(-self.expr, **self.assumptions)

    def _incompatible(self, x, op, reason=''):

        raise ValueError("Cannot determine %s(%s) %s %s(%s)%s" %
                         (self.__class__.__name__, self, op,
                          x.__class__.__name__, x, reason))

    def _incompatible_domains(self, x, op):

        self._incompatible(x, op, ' since the domains are incompatible')

    def _incompatible_quantities(self, x, op):

        self._incompatible(x, op,
                           """ since the units of the result are unsupported.
As a workaround use x.as_expr() %s y.as_expr()""" % op)

    def _dubious_operation(self, x, op):

        if op == '/':
            operation = 'deconvolution'
        else:
            operation = 'convolution'

        warn('Dubious operation: %s %s %s; you probably should be using %s' %
             (self, op, x, operation))

    def _add_compatible_domains(self, x):

        return self.domain == x.domain

    def _mul_domain(self, x):

        return self.domain

    def _mul_compatible_domains(self, x):

        if self.domain == x.domain:
            return True

        if (self.is_constant_domain or x.is_constant_domain):
            return True

        return False

    def _mul_dubious(self, x):

        return False

    def _div_domain(self, x):

        return self.domain

    def _div_compatible_domains(self, x):

        if self.domain == x.domain:
            return True

        if (self.is_constant_domain or x.is_constant_domain):
            return True

        return False

    def _div_dubious(self, x):

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
                self._incompatible(
                    x, op, ' since the units %s are incompatible with %s' % (self.units, x.units))

        cls = self.__class__
        xcls = x.__class__

        if x.is_constant_domain and x.quantity == 'undefined':
            if state.loose_units or x.expr == 0:
                # Allow voltage(1) + 2 etc.
                return cls, self, x, assumptions
            if self.is_transfer:
                # Allow transfer(1) == 1
                return cls, self, x, assumptions

        if (isinstance(self, (LaplaceDomainExpression, ZDomainExpression)) or
                isinstance(x, (LaplaceDomainExpression, ZDomainExpression))):
            assumptions = self.assumptions.add(x)

        if self.is_constant_domain and self.quantity == 'undefined':
            return xcls, self, x, assumptions

        if self.quantity == x.quantity:
            if self.is_constant_domain:
                return xcls, self, x, assumptions
            if x.is_constant_domain:
                return cls, self, x, assumptions
            if self._add_compatible_domains(x):
                return cls, self, x, assumptions

        # For phasor comparisons...
        # FIXME, remove phasor stuff?
        if self.is_phasor_ratio_domain and x.is_angular_fourier_domain:
            return cls, self, cls(x), assumptions
        elif self.is_angular_fourier_domain and x.is_phasor_ratio_domain:
            return xcls, cls(self), x, assumptions
        elif self.is_angular_frequency_response_domain and x.is_angular_fourier_domain:
            return cls, self, cls(x), assumptions
        elif self.is_angular_fourier_domain and x.is_angular_frequency_response_domain:
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
            if isinstance(x, (tuple, list, dict)):
                raise ValueError(
                    'Cannot multiply %s by a tuple, list, or dict %s' % (self, x))
            x = expr(x)

        # The following only applies to the generic expression classes.
        if (self.__class__ == FourierDomainExpression and
                x.__class__ == TimeDomainExpression):
            validate(state.f_times_t,
                     'Converting %s times %s to t-domain' % (self, x))
            return TimeDomainExpression(self.expr * x.expr)
        elif (x.__class__ == FourierDomainExpression and
              self.__class__ == TimeDomainExpression):
            validate(state.t_times_f,
                     'Converting %s times %s to t-domain' % (self, x))
            return TimeDomainExpression(self.expr * x.expr)
        elif (self.__class__ == AngularFourierDomainExpression and
                x.__class__ == TimeDomainExpression):
            validate(state.w_times_t,
                     'Converting %s times %s to t-domain' % (self, x))
            return TimeDomainExpression(self.expr * x.expr)
        elif (x.__class__ == AngularFourierDomainExpression and
              self.__class__ == TimeDomainExpression):
            validate(state.t_times_w,
                     'Converting %s times %s to t-domain' % (self, x))
            return TimeDomainExpression(self.expr * x.expr)
        elif (self.__class__ == LaplaceDomainExpression and
                x.__class__ == TimeDomainExpression):
            validate(state.s_times_t,
                     'Converting %s times %s to t-domain' % (self, x))
            return TimeDomainExpression(self.expr * x.expr)
        elif (x.__class__ == LaplaceDomainExpression and
              self.__class__ == TimeDomainExpression):
            validate(state.t_times_s,
                     'Converting %s times %s to t-domain' % (self, x))
            return TimeDomainExpression(self.expr * x.expr)

        # Try to convert immittance to a constant so that can handle I(t) * Z
        if x.is_immittance:
            try:
                xunits = x.units
                x = x.as_constant()
                x.units = xunits
            except:
                pass

        if not self._mul_compatible_domains(x):
            self._incompatible_domains(x, '*')

        if self._mul_dubious(x):
            self._dubious_operation(x, '*')

        if self.is_transform_domain:
            assumptions = self.assumptions.convolve(self, x)
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

        domain = self._mul_domain(x)
        if domain is None:
            self._incompatible(x, '*')

        if self.is_constant_domain:
            cls = x._class_by_quantity(quantity)
        else:
            cls = self._class_by_quantity(quantity, domain)

        value = self.expr * x.expr
        result = cls(value, **assumptions)
        result.units = self.units * x.units
        return result

    def __rmul__(self, x):
        """Reverse multiply."""

        return self.__mul__(x)

    def __truediv__(self, x, floor=False):
        """True divide."""

        if not isinstance(x, Expr):
            x = expr(x)

        # Try to convert immittance to a constant so that can handle V(t) / Z
        if x.is_immittance:
            try:
                x = x.as_constant()
            except:
                pass

        if not self._div_compatible_domains(x):
            self._incompatible_domains(x, '/')

        if self._div_dubious(x):
            self._dubious_operation(x, '/')

        assumptions = self.assumptions.convolve(self, x)

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

        domain = self._div_domain(x)
        if domain is None:
            self._incompatible(x, '/')

        if self.is_constant_domain:
            cls = x._class_by_quantity(quantity)
        else:
            cls = self._class_by_quantity(quantity, domain)

        if floor:
            value = self.expr // x.expr
        else:
            value = self.expr / x.expr
        result = cls(value, **assumptions)
        result.units = self.units / x.units

        return result

    def __rtruediv__(self, x, floor=False):
        """Reverse true divide."""

        from .matrix import Matrix

        if isinstance(x, Matrix):
            if floor:
                return x // self.expr
            else:
                return x / self.expr

        if not isinstance(x, Expr):
            x = expr(x)

        return x.__truediv__(self, floor)

    def __floordiv__(self, x):
        """Floor divide."""

        return self.__truediv__(x, floor=True)

    def __rfloordiv__(self, x):
        """Floor divide."""

        return self.__rtruediv__(x, floor=True)

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
            raise ValueError('Cannot compute %s(%s) ** %s' %
                             (self.__class__.__name__, self, x))

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

        # This does not speed up the comparison.
        # if self.expr == x.expr:
        #    return True

        if isinstance(self.expr, sym.Eq):
            return sym.simplify(self.expr.lhs - x.expr.lhs) == 0 and \
                sym.simplify(self.expr.rhs - x.expr.rhs) == 0

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

    def _cached_laplace(self):

        # The _laplace attribute stores a cached Laplace transform
        # of the expression.  This object can then store the cached
        # poles and zeros etc.
        try:
            return self._laplace
        except AttributeError:
            self._laplace = self.laplace()
            return self._laplace

    def cancel_terms(self):
        """Simplify terms in expression individually by converting
        each to rational functions."""

        result = 0
        for term in self.expr.as_ordered_terms():
            result += sym.cancel(term)
        return self.__class__(result, **self.assumptions)

    def convolve(self, x, commutate=False, **assumptions):
        """Convolve self with x.

        y(t) = int_{taumin}^{taumax} self(tau) x(t - tau) d tau

        If `commutate` is True, swap order of functions in integral.

        The result is an unevaluated integral.  It can be evaluated using
        the `doit()` method.

        Note, this method not simplify the convolution integral if one
        of the functions contains a Dirac delta.  This can be done
        calling the `simplify_dirac_delta()` method followed by the
        `simplify()` method.

        """

        if self.domain != x.domain:
            self._incompatible_domains(x, 'convolve')

        x = expr(x)
        f1 = self.expr
        f2 = x.expr
        if commutate:
            f1, f2 = f2, f1

        dummyvar = tausym if self.is_time_domain else nusym

        taumin = -oo
        taumax = oo

        if commutate:
            if x.is_causal:
                taumax = self.var
            if self.is_causal:
                taumin = 0
        else:
            if x.is_causal:
                taumin = 0
            if self.is_causal:
                taumax = self.var

        if x.is_causal and self.is_causal:
            assumptions['causal'] = True

        result = sym.Integral(f1.subs(self.var, self.var - dummyvar) *
                              f2.subs(self.var, dummyvar),
                              (dummyvar, taumin, taumax))
        ret = self.__class__(result, **assumptions)
        ret.units = self.units * x.units * self.domain_units
        return ret

    def parallel(self, x):
        """Parallel combination."""

        cls, self, x, assumptions = self.__compat_add__(x, '|')
        x = cls(x)

        return cls(self.expr * x.expr / (self.expr + x.expr), **assumptions)

    def copy(self):
        """Copy the expression."""
        return self.__class__(self.expr, **self.assumptions)

    @property
    def conj(self):
        """Return complex conjugate."""

        return self.__class__(sym.conjugate(self.expr), **self.assumptions)

    def conjugate(self):
        """Return complex conjugate."""

        return self.__class__(sym.conjugate(self.expr), **self.assumptions)

    @property
    def real(self):
        """Return real part.

         Note, SymPy does not always extract the real part.  For example,
        `exp(t*(-1 - sqrt(5)*j))/(20*sqrt(5) - 20*j)` or
        `exp(j * f) / (1 - exp(j * f))`

        A work-around is to use `rationalize_denominator()` or
        `expand(complex=True)` first.

        """

        assumptions = self.assumptions.copy()
        assumptions['real'] = True

        expr = self.expr
        # This can make operations such as abs really slow.
        # Without it, we will sometimes get Re functions.
        # expr = expr.expand(complex=True)

        dst = self.__class__(symsimplify(sym.re(expr)), **assumptions)
        dst.part = 'real'
        return dst

    @property
    def re(self):
        """Return real part, see `real`"""
        return self.real

    @property
    def imag(self):
        """Return imaginary part.

        Note, SymPy does not always extract the imaginary part.  For example,
        `exp(t*(-1 - sqrt(5)*j))/(20*sqrt(5) - 20*j)` or
        `exp(j * f) / (1 - exp(j * f))`

        A work-around is to use `rationalize_denominator()` or
        `expand(complex=True)` first.
        """

        assumptions = self.assumptions.copy()
        if self.is_real:
            dst = self.__class__(0, **assumptions)
            dst.part = 'imaginary'
            return dst

        assumptions['real'] = True

        expr = self.expr
        # This can make operations such as abs really slow.
        # Without it, we will sometimes get Im functions.
        # expr = expr.expand(complex=True)

        dst = self.__class__(symsimplify(sym.im(expr)), **assumptions)
        dst.part = 'imaginary'
        return dst

    @property
    def im(self):
        """Return imaginary part, see `im`"""
        return self.imag

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

        if (self.var is None or self.has(sym.Derivative) or
                self.has(sym.Integral)):
            self.__ratfun = None
        else:

            try:
                # Note, this handles expressions that are products of
                # rational functions and arbitrary delays.
                self.__ratfun = Ratfun(self.expr, self.var)
            except:
                self.__ratfun = None
        return self.__ratfun

    @property
    def ba(self):
        """Return lists of numerator and denominator coefficients."""

        a = self.D.coeffs()
        b = self.N.coeffs()

        a0 = a[0]
        if a0 != 1:
            a = ExprList([ax / a0 for ax in a])
            b = ExprList([bx / a0 for bx in b])
        return b, a

    @property
    def a(self):
        """Return list of denominator coefficients."""

        b, a = self.ba
        return a

    @property
    def b(self):
        """Return list of numerator coefficients."""

        b, a = self.ba
        return b

    @property
    def K(self):
        """Return gain.  Note, this removes the units."""

        # TODO, fix units
        return self.N.coeffs()[0] / self.D.coeffs()[0]

    @property
    def N(self):
        """Return numerator of rational function.
        The denominator is chosen so that it is a polynomial.
        Note, this removes the units."""

        # TODO, fix units
        return self.numerator

    @property
    def D(self):
        """Return denominator of rational function.
        The denominator is chosen so that it is a polynomial.
        Note, this removes the units."""

        # TODO, fix units
        return self.denominator

    @property
    def numerator(self):
        """Return numerator of rational function.
        The denominator is chosen so that it is a polynomial.
        Note, this removes the units."""

        N, D = self.as_N_D()
        return N

    @property
    def denominator(self):
        """Return denominator of rational function.
        The denominator is chosen so that it is a polynomial.
        Note, this removes the units."""

        N, D = self.as_N_D()
        return D

    def as_numer_denom(self, monic_denominator=False, use_sympy=False):

        return self.as_N_D(monic_denominator, use_sympy)

    def rationalize_denominator(self):
        """Rationalize denominator by multiplying numerator and denominator by
        complex conjugate of denominator."""

        N = self.N
        D = self.D
        Dconj = D.conj
        Nnew = (N * Dconj).simplify()
        # Dnew = (D * Dconj).simplify()
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

        c, r = factor_const(self.sympy, self.var)
        return ConstantDomainExpression(c), self.__class__(r, **self.assumptions)

    def term_const(self):

        from .utils import term_const

        c, r = term_const(self.sympy, self.var)
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

        from sympy import DiracDelta

        if self.has(DiracDelta):
            # Convert abs(X(f) + DiracDelta(f)) to abs(X(f)) + DiracDelta(f)
            rest, deltas = self.separate_dirac_delta()
            if not rest.has(DiracDelta):
                return rest.magnitude + deltas
            warn('Magnitude of expression with Dirac delta may be invalid: ', self)

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
        """Return magnitude in dB.  If squared voltage, squared current,
        or power, this uses 10 * log10 otherwise 20 * log10."""

        from sympy import DiracDelta

        if self.has(DiracDelta):
            rest, deltas = self.separate_dirac_delta()
            if not rest.has(DiracDelta) and rest != 0:
                return rest.dB + deltas

        # Need to clip for a desired dynamic range?
        # Assume reference is 1.
        if self.is_power or self.is_squared:
            dst = 10 * log10(self.magnitude)
        else:
            dst = 20 * log10(self.magnitude)
        dst.part = 'magnitude'
        dst.units = dB
        return dst

    @property
    def phase(self):
        """Return phase in radians."""

        if self.is_time_domain or self.is_discrete_time_domain:
            raise ValueError(
                'Cannot determine phase of time-domain expression %s' % self)

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
        """Returns True if expression does not have any free symbols (compare with `is_unchanging`)."""

        return self.expr.free_symbols == set()

    @property
    def is_realizable(self):
        """Return True if the signal's Laplace transform = 0
           in limit as s -> oo.  This is equivalent to the Laplace
           transform being strictly proper.

        Note, with this definition `cos(t)` is considered realizable."""

        return self._cached_laplace().is_strictly_proper

    @property
    def is_stable(self):
        """Return True if all the poles of the signal's Laplace transform
        are in the LH plane.

        Note, `exp(-alpha * t)` is considered stable since the
        unilateral Laplace transform ignores the region where `t < 0`.

        See also is_marginally_stable.
        """

        poles = self._cached_laplace().poles(aslist=True)
        for pole in poles:
            # If poles not constant, assume they are in the negative half-plane.
            if pole.is_constant and pole.real >= 0:
                return False
        return True

    @property
    def is_marginally_stable(self):
        """Return True if all the poles of the signal's Laplace transform
        are in the LH plane or on the boundary.

        Note, `exp(-alpha * t)` is considered stable since the
        unilateral Laplace transform ignores the region where `t < 0`.

        See also is_stable.
        """

        poles = self._cached_laplace().poles(aslist=True)
        for pole in poles:
            # If poles not constant, assume they are in the negative half-plane.
            if pole.is_constant and pole.real > 0:
                return False
        return True

    @property
    def is_unchanging(self):
        """Returns True if expression does not have a domain variable (compare with `is_constant`)."""

        if self.var is None:
            return True

        return self.var not in self.expr.free_symbols

    def evaluate(self, arg=None):
        """Evaluate expression at arg.  `arg` may be a scalar or a vector.
        The result is of type float or complex.

        If arg is iterable, a NumPy array is returned.

        There can be only one or fewer undefined variables in the expression.
        This is replaced by `arg` and then evaluated to obtain a result.
        """

        import numpy as np

        is_time = self.is_time_domain or self.is_discrete_time_domain
        is_causal = is_time and self.is_causal

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
                    arg = 500

                return np.exp(arg)

            def rect(arg):
                # Define in terms of Heaviside for consistency
                return heaviside(arg + 0.5) - heaviside(arg - 0.5)

            def sign(arg):
                # Define in terms of Heaviside for consistency
                return 2 * heaviside(arg) - 1

            def dtsign(arg):
                # Define in terms of unitstep for consistency
                return 2 * unitstep(arg) - 1

            def dtrect(arg):
                # Define in terms of UnitStep for consistency
                return unitstep(arg + 0.5) - unitstep(arg - 0.5)

            def sinc(arg):
                """SymPy sinc."""

                # This is used for sinc created by sympify, e.g., a =
                # expr('sinc(n)').  SymPy uses the unnormalized form
                # but NumPy (and Lcapy) use the normalized form.
                # Lambdify does some jiggery pokery and divides the
                # arg by pi since it is expecting that the NumPy sinc
                # function is going to be used.  SymPy's choice is
                # unfortunate from a numerical accuracy point of view
                # since sincn(n) should be zero for integer n, n != 0.

                # Undo SymPy jiggery pokery.
                arg = arg * np.pi

                return 1.0 if arg == 0 else np.sin(np.pi * arg) / (np.pi * arg)

            def sincn(arg):
                """Normalized sinc."""

                # Note, if sincn is made to print sinc, then lambdify will
                # call sinc.   Grrrr.

                return 1.0 if arg == 0 else np.sin(np.pi * arg) / (np.pi * arg)

            def sincu(arg):
                """Unnormalized sinc."""

                return 1.0 if arg == 0 else np.sin(arg) / arg

            def psinc(M, arg):
                """Periodic sinc."""

                D = np.sin(np.pi * arg)
                return 1.0 if D == 0 else np.sin(M * np.pi * arg) / (M * D)

            def trap(arg, alpha):

                absarg = abs(arg)
                foo = absarg - 0.5

                if alpha == 0:
                    if arg < -0.5 or arg > 0.5:
                        return 0.0
                    return 1.0

                if foo >= 0.5 * alpha:
                    return 0.0
                elif foo <= -0.5 * alpha:
                    return 1.0
                else:
                    return 0.5 - foo / alpha

            def tri(arg):

                if arg >= 1:
                    return 0.0
                elif arg <= -1:
                    return 0.0
                else:
                    return 1.0 - abs(arg)

            def ramp(arg):

                if arg >= 0:
                    return arg
                return 0

            def rampstep(arg):

                if arg < 0:
                    return 0
                elif arg > 1:
                    return 1
                else:
                    return arg

            def dirac(arg):
                return np.inf if arg == 0.0 else 0.0

            def unitimpulse(arg):
                return 1.0 if arg == 0 else 0.0

            def unitstep(arg, zero=None):
                if arg == 0:
                    if zero is None:
                        zero = unitstep_zero
                    return zero
                return 1.0 if arg >= 0 else 0.0

            def heaviside(arg, zero=None):
                if arg == 0:
                    if zero is None:
                        zero = heaviside_zero
                    return zero
                return 1.0 if arg > 0.0 else 0.0

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
                             [{'DiracDelta': dirac,
                              'Heaviside': heaviside,
                               'UnitImpulse': unitimpulse,
                               'UnitStep': unitstep,
                               'dtrect': dtrect, 'dtsign': dtsign,
                               'sinc': sinc, 'sincn': sincn,
                               'sincu': sincu, 'psinc': psinc,
                               'rect': rect, 'tri': tri, 'trap': trap,
                               'ramp': ramp, 'rampstep': rampstep,
                               'sqrt': sqrt, 'exp': exp, 'sign': sign},
                              "scipy", "numpy", "math", "sympy"])

            def func(arg):
                # Lambdify barfs on (-1)**n if for negative values of n.
                # even if have (-1)**n * Heaviside(n)
                # So this function heads Lambdify off at the pass,
                # if the function is causal.

                if is_causal and arg < 0:
                    return 0
                try:
                    result = func1(arg)
                except ZeroDivisionError:
                    result = complex(expr.limit(var, arg))

                # If get NaN evaluate limit.  This helps for sin(t) / t.
                if np.isnan(result):
                    result = complex(expr.limit(var, arg))
                # u(t) -
                if np.isinf(result):
                    result = complex(sym.simplify(expr).limit(var, arg))
                return result

            try:
                # Try to flush out weirdness using first argument
                response = func(arg0)
            except NameError as e:
                raise RuntimeError(
                    'Cannot evaluate expression %s: %s' % (self, e))
            except AttributeError as e:
                # Could return NaN but this would take some jiggery
                # pokery.  One solution is to find Piecewise with a
                # single clause, such as t >= 0 add another Piecewise
                # clause for t < 0 to return NaN.  Adding Piecewise to
                # lambdify functions does not seem to work.
                if expr.is_Piecewise:
                    raise RuntimeError(
                        'Cannot evaluate expression %s,'
                        ' due to undetermined conditional result' % self)

                raise RuntimeError(
                    'Cannot evaluate expression %s,'
                    ' probably have a mysterious function: %s' % (self, e))

            except TypeError as e:
                raise RuntimeError(
                    'Cannot evaluate expression %s: %s' % (self, e))

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

        # Use doit to expand Sum, etc.
        expr = self.doit().expr

        if not hasattr(self, 'var') or self.var is None:
            symbols = list(expr.free_symbols)
            if arg is None:
                if len(symbols) == 0:
                    return expr.evalf()
                raise ValueError(
                    'Undefined symbols %s in expression %s' % (tuple(symbols), self))
            if len(symbols) == 0:
                print('Ignoring arg %s' % arg)
                return expr.evalf()
            elif len(symbols) == 1:
                return evaluate_expr(expr, symbols[0], arg)
            else:
                raise ValueError(
                    'Undefined symbols %s in expression %s' % (tuple(symbols), self))

        var = self.var
        # Use symbol names to avoid problems with symbols of the same
        # name with different assumptions.
        varname = var.name
        free_symbols = set([symbol.name for symbol in expr.free_symbols])
        if varname in free_symbols:
            free_symbols -= set((varname, ))
            if free_symbols != set():
                raise ValueError('Undefined symbols %s in expression %s' % (
                    tuple(free_symbols), self))

        if arg is None:
            if expr.has(var):
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

        tweak_patterns = [pattern.expr if isinstance(
            pattern, (Expr, Function)) else pattern for pattern in patterns]
        return self.expr.has(*tweak_patterns)

    def has_symbol(self, sym):
        """Test if have symbol contained.  For example,
        V.has_symbol('a')
        V.has_symbol(t)

        """

        return self.has(expr(sym))

    def _subs1(self, old, new, domain=None, **kwargs):

        # This will fail if a variable has different attributes,
        # such as positive or real.

        # Should check for bogus substitutions, such as t for s.
        # Probably should disallow domain changing. Will need to
        # iterate through the domain variables and check if new
        # contains one.  Then disallow unless same as self.var.  Will
        # also need to handle s-> jw in select() and transform().

        # Allow domain=s as well as domain=`laplace`
        if domain is not None and not isinstance(domain, str):
            domain = domain.domain

        if new is old:
            return self

        expr = new
        if isinstance(new, Expr):
            if old == self.var or self.var is None:
                if domain is None:
                    domain = new.domain
            else:
                if domain is None:
                    domain = self.domain
            expr = new.expr
        else:
            domain = self.domain
            expr = sympify(expr)

        old = symbol_map(old)
        if isinstance(old, Expr):
            old = old.expr

        if isinstance(expr, list):
            # Get lists from solve.  These stymie sympy's subs.
            if len(expr) == 1:
                expr = expr[0]
            else:
                warn('Substituting a list...')

        if isinstance(old, str):
            result = self.sympy
        else:
            if state.warn_subs and not self.sympy.has(old):
                warn('Expression %s does not have %s' % (self.sympy, old))
            if state.break_subs and not self.sympy.has(old):
                import pdb
                pdb.set_trace()

            result = self.sympy.subs(old, expr, **kwargs)

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
        from numpy import ndarray, array

        if isinstance(arg, (tuple, list)):
            return [self._subs1(self.var, arg1) for arg1 in arg]

        if isinstance(arg, ndarray):
            return array([self._subs1(self.var, arg1) for arg1 in arg])

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

    def separate_dirac_delta(self):
        """Separate Dirac delta terms from expression."""

        from sympy import DiracDelta
        from .utils import separate_dirac_delta

        rest, deltas = separate_dirac_delta(self.sympy)
        rest = expr(rest)

        if rest.has(DiracDelta):
            rest, deltas = separate_dirac_delta(self.sympy.expand())

        rest = expr(rest)
        if rest.has(DiracDelta):
            warn('Could not separate all the Dirac deltas: %s' % self)

        deltaexpr = 0
        for delta in deltas:
            deltaexpr += delta

        return self.__class__(rest), self.__class__(deltaexpr)

    def simplify(self, **kwargs):
        """Simplify expression.

        This throws the kitchen sink at the problem but can be slow.

        See also simplify_terms and simplify_factors."""

        if self.has(AppliedUndef) and not \
           (self.has(sym.Integral) or self.has(sym.Derivative)):
            # This might be dodgy...  It replaces v(t) with v to help
            # simplification and then replaces v with v(t) at the end.
            # It cannot be used for Integral or Derivative since
            # the dependence is lost.

            new, defs = self.remove_undefs(return_mappings=True)
            return new.simplify(**kwargs).subs(defs)

        ret = symsimplify(self.expr, **kwargs)
        return self.__class__(ret, **self.assumptions)

    def simplify_units(self):
        """Simplify units into canonical form."""

        ret = self.__class__(self, **self.assumptions)
        ret.units = units.simplify_units(self.units)
        return ret

    def simplify_conjugates(self, **kwargs):
        """Combine complex conjugate terms."""

        result = simplify_conjugates(self.expr)
        return self.__class__(result, **self.assumptions)

    def simplify_factors(self, **kwargs):
        """Simplify factors in expression individually."""

        result = 0
        for factor in self.expr.as_ordered_factors():
            result *= symsimplify(factor, **kwargs)
        return self.__class__(result, **self.assumptions)

    def simplify_terms(self, **kwargs):
        """Simplify terms in expression individually."""

        result = 0
        for term in self.expr.as_ordered_terms():
            result += symsimplify(term, **kwargs)
        return self.__class__(result, **self.assumptions)

    def simplify_sin_cos(self, as_cos=False, as_sin=False):
        """Simplify c * cos(theta) - s * sin(theta) as A * cos(theta - phi)."""

        result = simplify_sin_cos(self.expr, as_cos, as_sin)
        return self.__class__(result, **self.assumptions)

    def simplify_dirac_delta(self):
        """Simplify DiracDelta(4 * t + 2) to DiracDelta(t + 0.5) / 4
        and DiracDelta(t) * x(t) to DiracDelta * x(0)."""

        result = simplify_dirac_delta(self.expr, self.var)
        return self.__class__(result, **self.assumptions)

    def simplify_heaviside(self):
        """Simplify Heaviside(4 * t + 2) to Heaviside(t + 0.5)
        and Heaviside(t)**2 to Heaviside(t), etc."""

        result = simplify_heaviside(self.expr, self.var)
        return self.__class__(result, **self.assumptions)

    def simplify_unit_impulse(self):
        """Simplify UnitImpulse(4 * k + 8) to UnitImpulse(k + 2), etc."""

        result = simplify_unit_impulse(self.expr, self.var)
        return self.__class__(result, **self.assumptions)

    def simplify_rect(self):
        """Simplify rect(4 * t + 2) to rect(t + 0.5)
        and rect(t)**2 to rect(t), etc."""

        result = simplify_rect(self.expr, self.var)
        return self.__class__(result, **self.assumptions)

    def expand_hyperbolic_trig(self):
        """Convert cosh(x) to exp(x) + exp(-x), etc."""

        result = expand_hyperbolic_trig(self.expr)
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
        """Substitute symbols in expression.

        `args` is either:
        - two arguments, e.g. foo.subs(old, new)
        - one iterable argument, e.g. foo.subs(iterable). The iterable may be
           o an iterable container with (old, new) pairs. In this case the
             replacements are processed in the order given with successive
             patterns possibly affecting replacements already made.
           o a dict or set whose key/value items correspond to old/new pairs.

        The domain of the result can be coerced using `domain` argument."""

        if len(args) > 2:
            raise ValueError('Too many arguments')
        if len(args) == 0:
            raise ValueError('No arguments')

        if len(args) == 2:
            return self._subs1(args[0], args[1], **kwargs)

        if isinstance(args[0], (dict, set)):
            # TODO, sort as per SymPy.
            dst = self
            for key, val in args[0].items():
                dst = dst._subs1(key, val, **kwargs)
            return dst
        elif isinstance(args[0], (list, tuple)):
            dst = self
            for key, val in args[0]:
                dst = dst._subs1(key, val, **kwargs)
            return dst

        if self.var is None:
            return self
        return self._subs1(self.var, args[0], **kwargs)

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
            if self.domain_units != 1:
                label += ' (%s)' % self.domain_units
        return label

    def differentiate(self, arg=None):
        """Differentiate expression."""

        if arg is None:
            arg = self.var

        arg = expr(arg)
        sarg = delcapify(arg)

        result = self.__class__(sym.diff(self.expr, sarg), **self.assumptions)
        result.units /= arg.units
        return result

    def diff(self, arg=None):

        return self.differentiate(arg)

    def doit(self, **hints):
        """Evaluate unevaluated functions such as integrals and sums."""

        result = self.__class__(self.expr.doit(**hints), **self.assumptions)
        result.part = self.part
        return result

    def integrate(self, arg=None, **kwargs):
        """Integrate expression.

        For example `exp(-3 * t).integrate((t, 0, oo))` gives `1 / 3`.

        """

        if arg is None:
            arg = self.var

        arg = expr(arg)
        sarg = delcapify(arg)

        result = self.__class__(sym.integrate(self.expr, sarg, **kwargs),
                                **self.assumptions)
        if isinstance(arg, tuple):
            arg = arg[0]

        result.units *= arg.units
        return result

    def rewrite(self, *args, **hints):
        """Rewrite expression.

        For example, `exp(j*a*t).rewrite(cos)` gives `ⅉ⋅sin(4⋅t) +
        cos(4⋅t)`.  Similarly, `cos(2 * t).rewrite(exp)` will expand
        the cosine as two complex exponentials."""

        args = delcapify(args)
        return self.__class__(self.sympy.rewrite(*args, **hints),
                              **self.assumptions)

    def solve(self, *symbols, **flags):
        """Solve expression.  This returns a list of solutions.  An empty list
        is returned if no solutions are found.  Note, by default,
        Lcapy assumes symbols are positive so an innocuous expression may fail
        to give a result if the solution is negative.

        For example:
        `>>> x = symbols('x')
         >>> y = x + 3
         >>> y.solve(x)
         []
         >>> x = symbols('x', positive=False)
         >>> y = x + 3
         >>> y.solve(x)
         [-3]`

        If no symbols are specified, all free symbols are solved for.

        See also `nsolve` for a numerical solver.
        """

        if self.has(AppliedUndef):
            new, defs = self.remove_undefs(return_mappings=True)
            return new.solve(*symbols, **flags).subs(defs)

        symbols = delcapify(ExprTuple(symbols).remove_undefs())
        symbols = [symbol_map(symbol) for symbol in symbols]
        solutions = sym.solve(self.expr, *symbols, **flags)
        return expr(solutions)

    def nsolve(self, x0=0, **kwargs):
        """Solve expression numerically.  This returns a list of solutions.
        An empty list is returned if no solutions are found.

        `x0` is starting point close to the solution.

        See sympy.solvers.solvers.nsolve for details.

        See also `solve` for a symbolic solver.

        """
        from sympy.solvers.solvers import nsolve

        return expr(nsolve(self.expr, x0, **kwargs))

    def split_dirac_delta(self):
        """Return expression as a list of terms.  The first term has no
        DiracDeltas, the second term collates the DiracDeltas, the
        third term collates derivatives of DiracDeltas, etc.

        For example, u(t) + DiractDelta(t, 1) returns
        [u(t), 0, DiracDelta(t, 1)]

        """

        # TODO: wrap return value as ExprList
        return split_dirac_delta(self)

    @property
    def symbols(self):
        """Return dictionary of symbols in the expression keyed by name."""

        expr = self.sympy
        symdict = {sym.name: sym for sym in expr.free_symbols}

        # Look for V(s), etc.
        funcdict = {
            atom.func.__name__: atom for atom in expr.atoms(AppliedUndef)}

        symdict.update(funcdict)
        return symdict

    def _fmt_roots(self, roots, aslist=False, pairs=False):

        units = uu.rad / uu.s

        def _wrap_dict(roots):

            rootsdict = {}
            for root, n in roots.items():
                rootsdict[expr(root, units=units)] = n
            return expr(rootsdict)

        def _wrap_list(roots):

            rootslist = []
            for root, n in roots.items():
                rootslist += [expr(root, units=units)] * n
            return expr(rootslist)

        if pairs:
            pairs, singles = pair_conjugates(roots)
            if aslist:
                return _wrap_list(pairs), _wrap_list(singles)
            else:
                return _wrap_dict(pairs), _wrap_dict(singles)

        if aslist:
            return _wrap_list(roots)
        else:
            return _wrap_dict(roots)

    def roots(self, aslist=False, pairs=False):
        """Return roots of expression as a dictionary.  Note this may not find
        them all.  In particular, if the rational function has a
        degree of five or higher.

        If `pairs` is True, return two dictionaries.  The first
        contains the conjugate pairs and the second contains the
        others

        If `aslist` is True, return roots as list.

        """

        if self._ratfun is None:
            roots = {}
        else:
            roots = self._ratfun.roots()
        return self._fmt_roots(roots, aslist, pairs)

    def zeros(self, aslist=False, pairs=False):
        """Return zeros of expression as a dictionary Note this may not find
        them all.  In particular, if the denominator of the rational
        function has a degree of five or higher.

        If `pairs` is True, return two dictionaries.  The first
        contains the conjugate pairs and the second contains the
        others

        If `aslist` is True, return zeros as list.

        """

        if self._ratfun is None:
            return self.N.roots(aslist, pairs)

        try:
            zeros = self._zeros
        except AttributeError:
            zeros = self._ratfun.zeros()
            self._zeros = zeros

        return self._fmt_roots(zeros, aslist, pairs)

    def poles(self, aslist=False, damping=None, pairs=False):
        """Return poles of expression as a dictionary.  Note this may not find
        them all.  In particular, if the denominator of the rational
        function has a degree of five or higher.

        If `pairs` is True, return two dictionaries.  The first
        contains the conjugate pairs and the second contains the
        others.

        If `aslist` is True, return poles as list.

        """

        if self._ratfun is None:
            # Handle expressions such as A(s) * (1 - exp(-s * T)) / B(s).
            return self.D.roots(aslist, pairs)

        try:
            poles = self._poles
        except AttributeError:
            poles = self._ratfun.poles(damping=damping)
            self._poles = poles

        polesdict = {}
        for pole in poles:
            key = pole.expr
            if key in polesdict:
                polesdict[key] += pole.n
            else:
                polesdict[key] = pole.n

        return self._fmt_roots(polesdict, aslist, pairs)

    def _as_ZPK(self):

        if self._ratfun is None:
            return None, None, None, None

        zeros, poles, K, undef = self._ratfun.as_ZPK()

        return zeros, poles, K, undef

    def parameterize_ZPK(self, zeta=None, ZPK=None):

        def def1(defs, symbolname, value, units):

            sym1 = symbol(symbolname, override=False)
            defs[symbolname] = expr(value, units=units)
            return sym1

        zeros, poles, K, undef = self._as_ZPK()
        if zeros is None and poles is None and K is None:
            return self, {}

        radpers = uu.rad / uu.s

        defs = ExprDict()
        K = def1(defs, 'K', K * undef, self.units *
                 radpers**(len(poles) - len(zeros)))

        N = 1
        D = 1
        for m, zero in enumerate(zeros):
            z = def1(defs, 'z%d' % (m + 1), zero, radpers)
            N *= sym.Add(self.var, -z.sympy, evaluate=False)

        for m, pole in enumerate(poles):
            p = def1(defs, 'p%d' % (m + 1), pole, radpers)
            D *= sym.Add(self.var, -p.sympy, evaluate=False)

        result = sym.Mul(K.sympy, sym.Mul(N, sym.Pow(D, -1)), evaluate=False)
        return self.__class__(result, **self.assumptions), defs

    def parameterize(self, zeta=None, ZPK=None):
        """Parameterize first and second-order expressions.

        For example, pexpr, defs = expr.parameterize()

        If parameterization is successful, defs is a dictionary
        of the parameters.  The original expression can be obtained
        with pexpr.subs(defs)

        For first order systems, parameterize as:

        K * (s + beta) / (s + alpha)

        K / (s + alpha)

        K (s + beta)

        where appropriate.

        If `zeta` is True, parameterize second-order expression in
        standard form using damping factor and natural frequency
        representation, i.e.

        N(s) / (s**2 + 2 * zeta * omega_0 * s + omega_0**2)

        otherwise parameterize as

        N(s) / (s**2 + 2 * sigma_1 * s + omega_1**2 + sigma_1**2)

        """

        def def1(defs, symbolname, value, units):
            # from .cexpr import cexpr

            sym1 = symbol(symbolname, override=False)
            defs[symbolname] = expr(value, units=units)
            return sym1

        if zeta is None and ZPK is None:
            zeta = True
            ZPK = False

        if ZPK:
            return self.parameterize_ZPK()

        factors = self.as_ordered_factors()

        spowers = [s**-4, s**-3, s**-2, s**-1, s, s**2, s**3, s**4]
        for spower in spowers:
            if spower in factors:
                result, defs = (self / spower).parameterize(zeta)
                return result * spower, defs

        radpers = uu.rad / uu.s
        N = self.N
        D = self.D

        ndegree = N.degree
        ddegree = D.degree
        ncoeffs = N.coeffs(norm=True)
        dcoeffs = D.coeffs(norm=True)

        result = None
        defs = ExprDict()

        K = self.K
        if ndegree < 1 and ddegree < 1:
            result = self
        elif ndegree == 1 and ddegree == 1:
            K = def1(defs, 'K', K, sym.S.One)
            alpha = def1(defs, 'alpha', dcoeffs[1], radpers)
            beta = def1(defs, 'beta', ncoeffs[1], radpers)
            result = K * (s + beta) / (s + alpha)
        elif ndegree == 1 and ddegree == 0:
            K = def1(defs, 'K', K, self.units / radpers)
            beta = def1(defs, 'beta', ncoeffs[1], radpers)
            result = K * (s + beta)
        elif ndegree == 0 and ddegree == 1:
            K = def1(defs, 'K', K, self.units * radpers)
            alpha = def1(defs, 'alpha', dcoeffs[1], radpers)
            result = K / (s + alpha)
        elif ddegree == 2:
            K = def1(defs, 'K', K, self.units * radpers**(ddegree - ndegree))
            coeffs = self.N.coeffs()

            if not zeta:
                sigma1 = def1(defs, 'sigma_1', dcoeffs[1] / 2, radpers)
                omega1 = def1(defs, 'omega_1',
                              sqrt(dcoeffs[2] - (dcoeffs[1] / 2)**2).simplify(), radpers)
                result = K * (self.N / coeffs[0]) / (s**2 + 2 *
                                                     sigma1 * s + sigma1**2 + omega1**2)
            else:
                omega0 = def1(defs, 'omega_0', sqrt(dcoeffs[2]), radpers)
                zeta = def1(defs, 'zeta',
                            dcoeffs[1] / (2 * sqrt(dcoeffs[2])), sym.S.One)
                result = K * (self.N / coeffs[0]) / \
                    (s**2 + 2 * zeta * omega0 * s + omega0**2)

        if result is None:
            # Copy?
            result = self

        return self.__class__(result, **self.assumptions), defs

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
        if self.is_Equality:
            return equation(self.lhs.canonical(), self.rhs.canonical())

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

        if self.is_Equality:
            return equation(self.lhs.general(), self.rhs.general())

        if self._ratfun is None:
            return self.copy()
        return self.__class__(self._ratfun.general(), **self.assumptions)

    def partfrac(self, combine_conjugates=False, pairs=False, damping=None,
                 method=None):
        """Convert rational function into partial fraction form.   For example,

        5 + (5 - 15 * j / 4) / (s + 2 * j) + (5 + 15 * j / 4) / (s - 2 * j)

        If `combine_conjugates` or `pairs` is True then the pair of partial
        fractions for complex conjugate poles are combined.   This creates
        a sum of biquad sections.

        `method` can be 'sub' (substitution method, the default) or
        'ec' (equating cofficients method).

        See also canonical, standard, general, timeconst, and ZPK."""

        pairs = pairs or combine_conjugates

        if self.is_Equality:
            return equation(self.lhs.partfrac(pairs, damping, method),
                            self.rhs.partfrac(pairs, damping, method))

        try:
            if self._ratfun is None:
                return self.copy()
            return self.__class__(self._ratfun.partfrac(pairs, damping, method),
                                  **self.assumptions)
        except ValueError:
            return self.as_sum().partfrac(pairs, damping, method)

    def recippartfrac(self, combine_conjugates=False, pairs=False,
                      damping=None, method=None):
        """Convert rational function into partial fraction form
        using reciprocal of variable.

        For example, if H = 5 * (s**2 + 1) / (s**2 + 5*s + 4)
        then H.recippartfrac() gives
        5/4 - 10/(3*(1 + 1/s)) + 85/(48*(1/4 + 1/s))

        If combine_conjugates is True then the pair of partial
        fractions for complex conjugate poles are combined.

        `method` can be 'sub' (substitution method, the default) or
        'ec' (equating cofficients method).

        See also canonical, standard, general, partfrac, timeconst, and ZPK."""

        pairs = pairs or combine_conjugates

        if self.is_Equality:
            return equation(self.lhs.recippartfrac(pairs, damping, method),
                            self.rhs.recippartfrac(pairs, damping, method))

        if self._ratfun is None:
            return self.copy()

        tmpsym = miscsymbol('qtmp')

        expr = self.subs(1 / tmpsym)
        ratfun = Ratfun(expr.expr, tmpsym)

        nexpr = ratfun.partfrac(combine_conjugates, damping, method)
        nexpr = nexpr.subs(tmpsym, 1 / self.var)

        return self.__class__(nexpr, **self.assumptions)

    def standard(self):
        """Convert rational function into mixed fraction form.  For example,

        (5 * s - 5) / (s**2 + 4) + 5

        This is the sum of strictly proper rational function and a
        polynomial.

        See also canonical, general, partfrac, timeconst, and ZPK.

        """

        if self.is_Equality:
            return equation(self.lhs.standard(), self.rhs.standard())

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

        if self.is_Equality:
            return equation(self.lhs.timeconst(), self.rhs.timeconst())

        if self._ratfun is None:
            return self.copy()
        return self.__class__(self._ratfun.timeconst(), **self.assumptions)

    def timeconst_terms(self):
        """Convert each term of expression into time constant form."""

        if self.is_Equality:
            return equation(self.lhs.timeconst_terms(), self.rhs.timeconst_terms())

        result = 0
        for term in self.expr.as_ordered_terms():
            result += self.__class__(term).timeconst()
        return self.__class__(result, **self.assumptions)

    def ZPK(self, pairs=False, combine_conjugates=False):
        """Convert to zero-pole-gain (ZPK) form (factored form).  For example,

        5 * (s + 1)**2 / ((s - 2 * j) * (s + 2 * j))

        If `combine_conjugates` or `pairs` is True, then conjugate pairs are
        combined to create a product of biquad sections.  For example,

        5 * (s + 1)**2 / (s**2 + 4)

        Note, both the numerator and denominator are expressed as
        products of monic factors, i.e., (s + 1 / 3) rather than (3 * s + 1).

        See also canonical, general, standard, partfrac, and timeconst.

        """
        if self.is_Equality:
            return equation(self.lhs.ZPK(), self.rhs.ZPK())

        if self._ratfun is None:
            return self.copy()
        return self.__class__(self._ratfun.ZPK(combine_conjugates or pairs),
                              **self.assumptions)

    def factored(self, pairs=False):
        """Convert to factored form.  For example,

        5 * (s + 1)**2 / ((s - 2 * j) * (s + 2 * j))

        If `pairs` is True, then conjugate pairs are combined.  For example,

        5 * (s + 1)**2 / (s**2 + 4)

        This is an alias for ZPK.  See also canonical, general,
        standard, partfrac, and timeconst.

        """

        if self.is_Equality:
            return equation(self.lhs.factored(), self.rhs.factored())

        if self._ratfun is None:
            return self.copy()
        return self.__class__(self._ratfun.ZPK(pairs), **self.assumptions)

    def expandcanonical(self):
        """Expand in terms for different powers with each term
        expressed in canonical form.  For example,

        s / (s**2 + 4) + 5 / (s**2 + 4)

        See also canonical, general, partfrac, timeconst, and ZPK."""

        if self.is_Equality:
            return equation(self.lhs.expandcanonoical(), self.rhs.expandcanonical())

        if self._ratfun is None:
            return self.copy()
        return self.__class__(self._ratfun.expandcanonical(), **self.assumptions)

    def expand_functions(self):
        """Expand functions in expression.

        For example `rect(t)` -> `u(t + 1 / 2) - u(t - 1 / 2)`
        """

        return self.__class__(expand_functions(self.sympy, self.var),
                              **self.assumptions)

    def expand_response(self):
        """Expand expression into a sum of Lcapy responses, where each
        response has the form B(var) * expp(-var * T) / A(var).
        This helps with inverse Laplace transforms.

        """

        N = self.N
        D = self.D
        Nterms = N.expand().as_ordered_terms()

        result = sym.S.Zero
        for Nterm in Nterms:
            result += (Nterm / D)
        return self.__class__(result, **self.assumptions)

    def coeffs(self, var=None, norm=False):
        """Return list of coeffs assuming the expr is a polynomial in terms of
        `var`.  If `var` is None, the default variable is used.
        The highest powers come first.

        This will fail for a rational function.  Instead use
        expr.N.coeffs or expr.D.coeffs for numerator or denominator
        respectively.

        If `norm` is True, normalize coefficients so highest power is 1.

        """

        if self._ratfun is None:
            return expr([self])

        if var is None:
            var = self.var

        try:
            z = sym.Poly(self.expr, var)
        except:
            raise ValueError(
                'Use .N or .D attribute to specify numerator or denominator of rational function')

        c = z.all_coeffs()
        if norm:
            return expr([sym.simplify(c1 / c[0]) for c1 in c])

        return expr(c)

    def normcoeffs(self, var=None):
        """Return list of coeffs (normalized so the highest power is 1)
        assuming the expr is a polynomial in s.  The highest powers
        come first.  This will fail for a rational function.  Instead
        use expr.N.normcoeffs or expr.D.normcoeffs for numerator or
        denominator respectively."""

        return self.coeffs(var, norm=True)

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

        result = self.copy()
        expr = self.expr
        result.expr = expr.replace(lambda expr: expr.is_Rational,
                                   lambda expr: sym.Float(expr))
        return result

    def floatrat(self):
        """This converts floating point numbers to rational numbers in an
        expression.  See also ratfloat.

        For example, 0.2 * t - > t / 5

        """

        expr = self.expr
        expr = expr.replace(lambda expr: expr.is_Float,
                            lambda expr: sym.sympify(str(expr), rational=True))

        return self.__class__(expr, **self.assumptions)

    def approximate_dominant(self, defs, threshold=0.01):
        """Approximate expression using ball-park numerical values for the
        symbols to decide which terms in a sum dominate the sum.  A term
        is neglected if its absolute value is below `threshold` times the
        maximum absolute value of the terms in the sum.

        `defs` is a dict, set, list, or tuple of symbol definitions;
        see `subs()`."""

        ndefs = {}
        for k, v in defs.items():
            ndefs[symbol_map(k)] = expr(v).sympy
        result = approximate_dominant(self.sympy, ndefs, threshold)
        return self.__class__(result, **self.assumptions)

    def approximate_exp(self, method='pade', order=1, numer_order=None):
        """Approximate exp(a).  The best time-domain response (without a jump)
        is achieved with 'numer_order == order - 1', where `numer` is
        the order of the numerator.  The best frequency-domain
        response is achieved with `numer_order == order`."""

        expr = approximate_exp(self, method, order, numer_order)
        return self.__class__(expr, **self.assumptions)

    def approximate_fractional_power(self, method='pade', order=2):
        """This is an experimental method to approximate
        s**a, where a is fractional, with a rational function using
        a Pade approximant of order `order`."""

        expr = approximate_fractional_power(self, method, order)
        return self.__class__(expr, **self.assumptions)

    def approximate_hyperbolic_trig(self, method='pade', order=1,
                                    numer_order=None):
        """Approximate cosh(a), sinh(a), tanh(a) with a rational function using
        a Pade approximant of order `order`."""

        expr = approximate_hyperbolic_trig(self, method, order, numer_order)
        return self.__class__(expr, **self.assumptions)

    def approximate_numer_order(self, order):
        """Reduce order of numerator to `order`.
        Warning, this assumes `self.var` is much smaller than 1."""

        N, D = self.as_N_D()

        expr = approximate_order(N.sympy, self.var, order) / D.sympy
        return self.__class__(expr, **self.assumptions)

    def approximate_denom_order(self, order):
        """Reduce order of denominator to `order`.
        Warning, this assumes `self.var` is much smaller than 1."""

        N, D = self.as_N_D()

        expr = N.sympy / approximate_order(D.sympy, self.var, order)
        return self.__class__(expr, **self.assumptions)

    def approximate_order(self, order):
        """Reduce order of numerator and denominator to `order`.
        Warning, this assumes `self.var` is much smaller than 1."""

        N, D = self.as_N_D()

        expr = approximate_order(N.sympy, self.var, order) / \
            approximate_order(D.sympy, self.var, order)
        return self.__class__(expr, **self.assumptions)

    def approximate(self, method='pade', order=1, numer_order=None, defs=None,
                    threshold=0.01):
        """Approximate an expression.  The order of attempted approximations is:
        1. Fractional powers
        2. Exponential functions
        3. Hyperbolic trigonometric functions
        4. Dominant terms (if `defs` is defined).

        See also `approximate_numer_order`, `approximate_denom_order`, and
        `approximate_order` methods.
        """

        result = self.approximate_fractional_power(method=method,
                                                   order=order)
        result = result.approximate_exp(method=method,
                                        order=order,
                                        numer_order=numer_order)
        result = result.approximate_hyperbolic_trig(method=method,
                                                    order=order,
                                                    numer_order=numer_order)
        if defs is not None:
            result = result.approximate_dominant(defs, threshold)
        return result

    def as_value_unit(self):
        """Return tuple of value and unit.  For example,

        >> > v = voltage(5)
        >> > v.as_value_unit
        (5, volts)
        """

        from .units import units

        return units.as_value_unit(self.expr)

    def as_N_D(self, monic_denominator=False, use_sympy=False):
        """Responses due to a sum of delayed transient responses
        cannot be factored into ZPK form with a constant delay.
        For example, sometimes SymPy gives:

            ⎛    s⋅τ     ⎞ - s⋅τ
            ⎝V₁⋅ℯ   -  V₂⎠⋅ℯ
        I = ────────────────────
               s⋅(L⋅s + R)

        This method tries to extract the numerator and denominator
        where the denominator is a polynomial.

        N, D = I.as_N_D()

                     -s⋅τ
        N = V₁ - V₂⋅ℯ
        D = s⋅(L⋅s + R)"""

        N, D = as_N_D(self.expr, self.var, monic_denominator, use_sympy)

        # Strip quantity and assumptions
        cls = self._class_by_quantity('undefined')
        return cls(N, **self.assumptions), cls(D, **self.assumptions)

    def as_sum(self):
        """Responses due to a sum of delayed transient responses
        cannot be factored into ZPK form with a constant delay.
        For example, sometimes SymPy gives:

            ⎛    s⋅τ     ⎞ - s⋅τ
            ⎝V₁⋅ℯ - V₂⎠⋅ℯ
        I = ────────────────────
               s⋅(L⋅s + R)

        While this cannot be factored into ZPK form, it can be
        expressed as a sum of ZPK forms or as a partial fraction
        expansion.  However, SymPy does not play ball if trying to
        express as a sum of terms:

        I.as_ordered_terms()

        ⎡⎛    s⋅τ     ⎞ - s⋅τ⎤
        ⎢⎝V₁⋅ℯ - V₂⎠⋅ℯ    ⎥
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
        """Convert expression into a continued fraction."""

        def foo(coeffs):

            if len(coeffs) == 1:
                return coeffs[0]
            return coeffs[0] + 1 / foo(coeffs[1:])

        coeffs = self.continued_fraction_coeffs()
        result = foo(coeffs)
        return self.__class__(result, **self.assumptions)

    def as_ratfun_delay(self):

        B, A, delay, undef = self._ratfun.as_B_A_delay_undef()
        if undef != 1:
            raise ValueError('Have undefined expression %s' % undef)

        return self.__class__(B / A, **self.assumptions), delay

    def _as_B_A_delay_undef(self):

        return self._ratfun.as_B_A_delay_undef()

    def continued_fraction_inverse_coeffs(self):
        """Convert expression into a continued fraction with inverse
        coefficients."""

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

    def estimate(self, x, y, method='trf', ranges=None, Ns=10, **kwargs):
        """Estimate the parameters of an expression that best fits measured data using
        non-linear least squares optimization techniques.

        `x` is an ndarray of values for the dependent variable,
        `y` is an ndarray of values for the independent variable,
        `method` is the optimization method (including: 'brute', 'trf' (default), 'dogbox',
        'Nelder-Mead`, `Powell`),
        `ranges` is a dictionary of the search ranges keyed by the parameter name,
        `Ns` is the number of steps per range for the 'brute' optimizer.

        The returned result is a `FitterResult` object.   This has attributes:
        `params` a dictionary of the estimated parameter values keyed by parameter name,
        `rmse` the root mean squared error.

        Here's an example of use:
        `
        >>> e = expr('a * exp(-t  / tau) * u(t)')
        >>> tv = arange(100)
        >>> vv = e.subs({'a': 1, 'tau': 10}).evaluate(tv)
        >>> results = e.estimate(tv, vv, ranges={'a': (0, 10), 'tau': (1, 20)})
        >>> results.params
        {'a': 1.000000000048109, 'tau': 9.999999998432187}
        >>> results.rmse
        3.489526384702217e-22
        `
        """

        from .fitter import fit

        return fit(self, x, y, method, ranges, **kwargs)

    def force_time(self):

        if state.force_time:
            return self.time()
        return self

    def remove_undefs(self, return_mappings=False):
        """Replace applied undefined functions(undefs) with symbols, for
        example, replace x(t) with x, etc.

        This is useful for simplifying and solving equations.

        This method gives up if it finds an expression such as
        x(n - 1) + 2 * x(n) since the arguments are different.

        If return_mappings is True, then a dictionary of substitutions
        is returned as well as the modified expression.  For example,

        new, defs = expr.remove_undefs(return_mappings=True)

        The original expression can be obtained using:

        new.subs(defs)

        """

        mappings = {}
        e = self.expr
        for item in sym.preorder_traversal(e):
            if isinstance(item, AppliedUndef):
                name = str(item)
                parts = name.split('(')
                name = parts[0]
                if name in mappings and mappings[name] != item:
                    # Have found something like x(n) + x(n - 1),
                    # so give up...   We could create different
                    # named symbols but these cannot conflict
                    # with other symbols in the expression.
                    break

                if name == 'I':
                    warn('Replacing %s with %s; this is interpreted as the imaginary number' % (
                        item, name))

                mappings[name] = item
                # Need to propagate complex assumption, etc.
                e = e.subs(item, expr(name).expr)

        ret = self.__class__(e, **self.assumptions)
        if return_mappings:
            return ret, mappings
        else:
            return ret

    def remove_condition(self):
        """Remove the piecewise condition from the expression.
        See also force_causal."""

        if not self.is_conditional:
            return self
        expr = self.expr
        expr = expr.args[0].args[0]
        return self.__class__(expr)

    def remove_images(self, m1=0, m2=0):
        """Remove all spectral images resulting from a DTFT.

        For example,

        >> > x = Sum(DiracDelta(f - m/Delta_t), (m, -oo, oo))
        >> > x.remove_images()
        DiracDelta(f)

        Alternatively, the number of images can be changed,
        for example,

        >> > x = Sum(DiracDelta(f - m/Delta_t), (m, -1, 1))
        >> > x.remove_images()
        Sum(DiracDelta(f - m/Delta_t), (f, -1, 1))

        """
        from .sym import dt

        var = self.var

        if var is fsym:
            scale = dt
        elif var is Fsym:
            scale = 1
        elif var is omegasym:
            scale = dt / (2 * pi)
        elif var is Omegasym:
            scale = 1 / (2 * pi)
        else:
            raise RuntimeError('Mystery var %s' % var)

        result = remove_images(self.expr, var, scale, m1, m2)
        return self.__class__(result, **self.assumptions)

    def as_QRD(self):

        return self._ratfun.as_QRD()

    def as_QRPO(self):

        return self._ratfun.as_QRPO()

    def zero_initial_conditions(self):

        if not self.has(AppliedUndef):
            return self

        result = 0
        expr = self.sympy
        terms = expr.as_ordered_terms()
        for term in terms:
            for factor in term.as_ordered_factors():
                if isinstance(factor, AppliedUndef) and factor.args[0] == 0:
                    term = 0
                    break
                elif isinstance(factor, sym.Subs) and factor.args[0].is_Derivative and \
                        isinstance(factor.args[0].args[0], AppliedUndef) and \
                        factor.args[2].args[0] == 0 and \
                        factor.args[1].args[0] == factor.args[0].args[0].args[0]:
                    term = 0
                    break

            result += term

        return self.__class__(result)


def exprcontainer(arg, **assumptions):

    from numpy import ndarray

    if isinstance(arg, (ExprList, ExprTuple, ExprDict)):
        return arg
    elif isinstance(arg, list):
        return ExprList(arg, **assumptions)
    elif isinstance(arg, (tuple, sym.Tuple)):
        return ExprTuple(arg, **assumptions)
    elif isinstance(arg, dict):
        return ExprDict(arg)
    elif isinstance(arg, ndarray):
        from .vector import Vector
        if arg.ndim > 1:
            raise ValueError(
                'Multidimensional arrays unsupported; convert to Matrix')
        return Vector(arg, **assumptions)

    raise ValueError('Unsupported exprcontainer %s' % arg.__class__.__name__)


def _make_domain(expr, **assumptions):

    symbols = expr.free_symbols

    if tsym in symbols:
        return texpr(expr, **assumptions)
    elif ssym in symbols:
        return sexpr(expr, **assumptions)
    elif fsym in symbols:
        return fexpr(expr, **assumptions)
    elif omegasym in symbols:
        return omegaexpr(expr, **assumptions)
    elif nsym in symbols:
        return nexpr(expr, **assumptions)
    elif ksym in symbols:
        return kexpr(expr, **assumptions)
    elif zsym in symbols:
        return zexpr(expr, **assumptions)
    elif Fsym in symbols:
        return Fexpr(expr, **assumptions)
    elif Omegasym in symbols:
        return Omegaexpr(expr, **assumptions)
    else:
        return cexpr(expr, **assumptions)


def expr(arg, override=False, units=None, **assumptions):
    """Create Lcapy expression from arg.

    If `arg` is an `Expr` it is returned, unless `assumptions` is specified.

    If `arg` is a string:
       If a t symbol is found in the string a TimeDomainExpression object is created.
       If a s symbol is found in the string a LaplaceDomainExpression object is created.
       If a f symbol is found in the string an FourierDomainExpression object is created.
       If an omega symbol is found in the string an AngularFourierDomainExpression object is created.

    For example, v = expr('3 * exp(-t / tau) * u(t)')

    V = expr('5 * s', causal=True)

    If `override` is True, then create new symbol(s) even if
    previously defined by SymPy.
    """

    # TODO: convert expr('x(t)') into Function('x')(t)

    from .sequence import Sequence

    if arg is None:
        return arg

    if isinstance(arg, Expr):
        if assumptions == {} and units is None:
            return arg
        new = arg.__class__(arg, **assumptions)
        if units is not None:
            new.units = units
        return new
    if isinstance(arg, Sequence):
        return arg

    if not isinstance(arg, str) and hasattr(arg, '__iter__'):
        return exprcontainer(arg, **assumptions)

    # Don't set rational=True since this will set rational
    # assumption for symbols.
    expr = sympify(arg, override=override, **assumptions)

    lexpr = _make_domain(expr, **assumptions)
    if not lexpr.has(uu.Quantity):
        if units is not None:
            lexpr.units = units
        return lexpr

    from .units import units

    cls = lexpr.__class__
    expr, units = units.as_value_unit(lexpr.expr)

    # 5 * t * u.volts -> V
    # 5 * cos(t) * u.volts -> V
    # 5 * s * u.volts -> V / Hz
    warn('This may be deprecated since the units may not be what you expect')

    if units == uu.volts:
        return cls(expr, **assumptions).as_voltage()
    elif units == uu.amperes:
        return cls(expr, **assumptions).as_current()
    elif units == uu.ohms:
        return cls(expr, **assumptions).as_impedance()
    elif units == uu.siemens:
        return cls(expr, **assumptions).as_admittance()
    elif units == uu.watts:
        return cls(expr, **assumptions).as_power()
    warn('Unhandled units: %s' % units)
    return lexpr


def expr_class(domain, arg, **assumptions):

    try:
        quantity = arg.quantity
    except:
        quantity = 'undefined'

    cls = expressionclasses.get_quantity(domain, quantity)
    return cls


def expr_make(domain, arg, **assumptions):

    cls = expr_class(domain, arg)
    return cls(arg, **assumptions)


def equation(lhs, rhs, inputsym='x', outputsym='y', **assumptions):
    """Create an Lcapy equation.

    This is an Lcapy expression of the form Eq(lhs, rhs).
    For example,
    e = equation('Y(s)', 'X(s) * 2 * s')

    The left hand side(lhs) and right hand side subexpressions
    can be obtained with the `lhs` and `rhs` attributes."""

    from .differenceequation import DifferenceEquation

    lhs = expr(lhs)
    rhs = expr(rhs)
    # Check if lhs and rhs compatible.
    diff = lhs - rhs

    if diff.is_discrete_time_domain:
        return DifferenceEquation(lhs, rhs, inputsym, outputsym, **assumptions)

    cls = lhs.__class__

    return cls(sym.Eq(lhs.expr, rhs.expr, evaluate=False), **assumptions)


def symbol(name, **assumptions):
    """Create Lcapy symbols from whitespace or comma delimited string of
    symbol names.

    By default, symbols are assumed to be positive and real.

    `symbols('a b', real=True)` creates real symbols

    `symbols('a b', complex=True)` creates complex symbols

    `symbols('a b', integer=True, positive=True)` creates positive integer symbols

    See also symbol.
    """

    ssym = usersymbol(name, **assumptions)
    # Create Lcapy symbol
    return expr(ssym, **assumptions)


def symbols(names, **assumptions):
    """Create Lcapy symbols from whitespace or comma delimited string of
    symbol names.  See also symbol."""

    from .parser import split

    namelist = split(names, ", ")
    symbols = []
    for name in namelist:
        symbols.append(symbol(name, **assumptions))
    if len(symbols) == 1:
        return symbols[0]
    return symbols


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


def delcapify(expr):
    """De-lcapify expression to create pure SymPy expression."""

    if isinstance(expr, tuple):
        return tuple([delcapify(arg) for arg in expr])
    elif isinstance(expr, list):
        return [delcapify(arg) for arg in expr]
    elif isinstance(expr, dict):
        ret = {}
        for key, val in expr.items():
            ret[delcapify(key)] = delcapify(val)
        return ret
    elif hasattr(expr, 'expr'):
        return expr.expr

    return expr


def check_expr(expr):

    args = getattr(expr, 'args', None)
    if args is not None:
        for arg in args:
            if isinstance(arg, Expr):
                print(arg)
            check_expr(arg)


from .cexpr import cexpr, ConstantDomainExpression  # nopep8
from .fexpr import f, fexpr, FourierDomainExpression  # nopep8
from .omegaexpr import omega, omegaexpr, AngularFourierDomainExpression  # nopep8
from .normfexpr import Fexpr  # nopep8
from .normomegaexpr import Omegaexpr  # nopep8
from .phasor import PhasorRatioDomainExpression  # nopep8
from .texpr import t, texpr, TimeDomainExpression  # nopep8
from .sexpr import s, sexpr, LaplaceDomainExpression  # nopep8
from .nexpr import nexpr  # nopep8
from .kexpr import kexpr  # nopep8
from .zexpr import zexpr, ZDomainExpression  # nopep8
from .expressionclasses import expressionclasses  # nopep8

# Horrible hack to work with IPython around Sympy's back for LaTeX
# formatting.  The problem is that Sympy does not check for the
# _repr_latex method and instead relies on a predefined list of known
# types.  See _can_print method in sympy/interactive/printing.py

try:
    from .printing import latex
    formatter = sys.displayhook.shell.display_formatter.formatters['text/latex']

    for cls in (ExprList, ExprTuple, ExprDict):
        formatter.type_printers[cls] = Expr._repr_latex_
except:
    pass
