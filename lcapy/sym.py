"""This module contains functions for handling symbols and in particular
for integrating with SymPy.

Copyright 2014--2022 Michael Hayes, UCECE

"""


from .config import excludes, aliases
from sympy.parsing.sympy_parser import parse_expr, auto_number, rationalize
try:
    from sympy.parsing.sympy_parser import NUMBER, NAME, OP
except:
    from sympy.parsing.sympy_tokenize import NUMBER, NAME, OP

from sympy import Basic, Symbol, Expr, Atom
from sympy.core.function import AppliedUndef
import sympy as sym
import re
from .state import state
from .simplify import simplify_dirac_delta, simplify_heaviside, simplify_rect
from .simplify import simplify_unit_impulse
from .extrafunctions import (UnitImpulse, UnitStep, rect, dtrect, tri,
                             trap, dtsign, sincn, sincu, psinc, ramp,
                             rampstep)
from warnings import warn

__all__ = ('sympify', 'simplify', 'symbol_delete')


global_dict = {}
exec('from sympy import *', global_dict)

for _alias, _name in aliases.items():
    global_dict[_alias] = global_dict[_name]

global_dict['abs'] = sym.Abs
# Populate Lcapy defined functions
global_dict['rect'] = rect
global_dict['tri'] = tri
global_dict['trap'] = trap
global_dict['sincn'] = sincn
global_dict['sincu'] = sincu
global_dict['psinc'] = psinc
global_dict['dtrect'] = dtrect
global_dict['dsign'] = dtsign
global_dict['ui'] = UnitImpulse
global_dict['us'] = UnitStep
global_dict['ramp'] = ramp
global_dict['rampstep'] = rampstep


for _symbol in excludes:
    try:
        global_dict.pop(_symbol)
    except:
        pass


symbol_dict = {}
for name in ['Symbol', 'Function', 'Integer']:
    symbol_dict[name] = global_dict[name]


def capitalize_name(name):

    return name[0].upper() + name[1:]


def symbol_name(symbol):

    return str(symbol)


def symbols_find(arg):
    """Return list of symbols in arg.  No symbols are cached."""

    symbols = []

    def find_symbol(tokens, local_dict, global_dict):

        for tok in tokens:
            tokNum, tokVal = tok
            if tokNum == NAME:
                name = tokVal

                # Hack to fix common naming problem.  Perhaps use a dictionary
                # of symbol aliases?
                if name == 'omega0':
                    name = 'omega_0'

                if name not in local_dict and name not in global_dict:
                    symbols.append(name)
        return ([(NUMBER, '0')])

    if isinstance(arg, str):
        parse_expr(arg, transformations=(find_symbol, ),
                   global_dict=global_dict, local_dict={}, evaluate=False)

        return symbols

    try:
        arg = arg.sympy
    except:
        pass

    if not isinstance(arg, (Symbol, Expr, AppliedUndef)):
        return []
    return [symbol_name(symbol) for symbol in arg.atoms(Symbol, AppliedUndef)]


def parse(string, symbols=None, evaluate=True, local_dict=None,
          global_dict=None, **assumptions):
    """Handle arbitrary strings that may refer to multiple symbols."""

    if symbols is None:
        symbols = {}

    if local_dict is None:
        local_dict = {}
    if global_dict is None:
        global_dict = {}

    cache = assumptions.pop('cache', True)

    def auto_symbol(tokens, local_dict, global_dict):
        """Inserts calls to ``Symbol`` or ``Function`` for undefined
        variables/functions."""
        result = []

        tokens.append((None, None))  # so zip traverses all tokens
        for tok, nextTok in zip(tokens, tokens[1:]):
            tokNum, tokVal = tok
            nextTokNum, nextTokVal = nextTok
            if tokNum == NAME:
                name = tokVal

                # Automatically add Function.  We ignore the assumptions.
                # These could be supported by modifying fourier.py/laplace.py
                # to propagate assumptions when converting V(s) to v(t), etc.
                if (nextTokVal == '(' and name not in global_dict):
                    result.extend([(NAME, 'Function'),
                                   (OP, '('), (NAME, repr(name)), (OP, ')')])
                    continue

                # Hack to fix common naming problem.  Perhaps use a dictionary
                # of symbol aliases?
                if name == 'omega0':
                    name = 'omega_0'

                if name in global_dict:

                    obj = global_dict[name]
                    if isinstance(obj, (Basic, type)):
                        result.append((NAME, name))
                        continue

                    if callable(obj):
                        result.append((NAME, name))
                        continue

                if name in local_dict:
                    # print('Found %s' % name)
                    # Could check assumptions.
                    result.append((NAME, name))
                    continue

                # Automatically add Symbol
                result.extend([(NAME, 'Symbol'),
                               (OP, '('), (NAME, repr(name))])
                for assumption, val in assumptions.items():
                    result.extend([(OP, ','),
                                   (NAME, '%s=%s' % (assumption, val))])
                result.extend([(OP, ')')])

            else:
                result.append((tokNum, tokVal))

        return result

    s = parse_expr(string, transformations=(auto_symbol, auto_number,
                                            rationalize),
                   global_dict=global_dict, local_dict=local_dict,
                   evaluate=evaluate)
    if not cache:
        return s

    # Add newly defined symbols.
    symbols.add_from_expr(s)
    return s


def sympify1(arg, symbols=None, evaluate=True, override=False, rational=True,
             **assumptions):
    """Create a SymPy expression.

    The purpose of this function is to head SymPy off at the pass and
    apply the defined assumptions.

    If `evaluate` is True, the expression is evaluated.

    If `override` is True, then create new symbol(s) even if
    previously defined by SymPy.
    """

    if symbols is None:
        symbols = {}

    if isinstance(arg, Expr):
        if not arg.has(sym.Float):
            return arg
        # This is needed to catch 0.1 + sym.I
        arg = arg.replace(lambda expr: expr.is_Float,
                          lambda expr: sym.sympify(str(expr), rational=rational))
        return arg

    if isinstance(arg, Symbol):
        return arg

    # Why doesn't SymPy do this?
    if isinstance(arg, complex):
        re = sym.sympify(str(arg.real), rational=rational, evaluate=evaluate)
        im = sym.sympify(str(arg.imag), rational=rational, evaluate=evaluate)
        if im == 1.0:
            arg = re + sym.I
        else:
            arg = re + sym.I * im
        return arg

    if isinstance(arg, float):
        # Note, need to convert to string to achieve a rational
        # representation.
        return sym.sympify(str(arg), rational=rational, evaluate=evaluate)

    if isinstance(arg, str):
        # Quickly handle the simple case.
        if arg in symbols:
            return symbols[arg]

        # Handle arbitrary strings that may refer to multiple symbols.
        if override:
            # Use restricted global symbol dictionary so that can
            # override pre-defined SymPy symbols.
            gdict = symbol_dict
        else:
            gdict = global_dict
        return parse(arg, symbols, evaluate=evaluate,
                     local_dict=symbols, global_dict=gdict, **assumptions)

    return sym.sympify(arg, rational=rational, locals=symbols,
                       evaluate=evaluate)


def sympify(expr, evaluate=True, override=False, rational=True, **assumptions):
    """Create a SymPy expression.

    By default, symbols are assumed to be positive reals if no assumptions
    are explicitly defined.

    Note, this will not modify previously defined symbols with the
    same name.  Thus you cannot change the assumptions.

    If `evaluate` is True, the expression is evaluated.

    If `override` is True, then create new symbol(s) even if
    previously defined by SymPy.
    """

    if 'positive' in assumptions:
        if not assumptions['positive']:
            assumptions.pop('positive')
    elif ('real' not in assumptions and 'complex' not in assumptions
          and 'integer' not in assumptions):
        assumptions['positive'] = True
        # Note this implies that imag is False.   Also note that all
        # reals are considered complex (but with a zero imag part).

    return sympify1(expr, state.symbols, evaluate, override,
                    rational, **assumptions)


def miscsymbol(name, force=False, **assumptions):
    """Create a SymPy symbol.

    This function allows symbol assumptions to be defined,
    e.g., `real=True, integer=True`

    By default, symbols are assumed to be positive unless real is
    defined.

    If `name` is already a symbol, it is overridden unless it is
    a domain symbol in which case `force` must be `True`.

    """

    return state.symbols.add_misc(name, force=force, **assumptions)


def usersymbol(name, force=False, **assumptions):
    """Create a SymPy symbol.

    This function allows symbol assumptions to be defined,
    e.g., `real=True, integer=True`

    By default, symbols are assumed to be positive unless real is
    defined.

    If `name` is already a symbol, it is overridden unless it is
    a domain symbol in which case `force` must be `True`.

    """

    try:
        name = name.sympy
    except:
        pass
    try:
        # Handle sympy Symbol
        name = name.name
    except:
        pass

    return state.symbols.add_user(name, force=force, **assumptions)


def domainsymbol(name: str, **assumptions):
    """Create a SymPy symbol and register as a domain symbol
    that should not be overwritten."""

    return state.symbols.add_domain(name, **assumptions)


def symsimplify(expr, var=None, **kwargs):
    """Simplify a SymPy expression.  This is a hack to work around
    problems with SymPy's simplify API."""

    # Handle Matrix types
    if hasattr(expr, 'applyfunc'):
        return expr.applyfunc(lambda x: symsimplify(x, **kwargs))

    if expr.has(sym.DiracDelta):
        expr = simplify_dirac_delta(expr)
    if expr.has(UnitImpulse):
        expr = simplify_unit_impulse(expr)
    if expr.has(sym.Heaviside, UnitStep):
        expr = simplify_heaviside(expr, var)
    if expr.has(rect, dtrect):
        expr = simplify_rect(expr, var)

    # This+ gets expanded into piecewise...
    if expr.has(sym.sign):
        # Could replace sign with dummy function, simplify, and then restore...
        # Could replace 1 + sign(t) with 2 * Heaviside(t) but this depends
        # on the definition of sign(t).  SymPy and NumPy define sign(0) = 0.
        # This implies Heaviside(0) = 0.5.
        return expr

    try:
        if expr.is_Function and expr.func in (sym.Heaviside, sym.DiracDelta):
            return expr
    except:
        pass

    expr = sym.simplify(expr, **kwargs)

    return expr


def simplify(expr, **kwargs):
    """Simplify an Lcapy expression.  This is not straightforward, see
    sympy.simplify."""

    try:
        return expr.simplify(**kwargs)
    except:
        pass

    try:
        expr = expr.sympy
    except:
        pass

    return symsimplify(expr, **kwargs)


def is_sympy(expr):
    return isinstance(expr, (Symbol, Expr, AppliedUndef))


def symdebug(expr, s='', indent=0):
    """See also the SymPy function srepr."""

    def _debug_args(args, s='', indent=0):

        for m, arg in enumerate(args):
            s = symdebug(arg, s, indent)
            if m == len(expr.args) - 1:
                s += ')\n'
            else:
                s += ',\n' + ' ' * indent
        return s

    if isinstance(expr, Symbol):
        s += str(expr) + ': %s' % expr.assumptions0

    elif isinstance(expr, Atom):
        s += str(expr)

    elif isinstance(expr, Expr):

        name = expr.__class__.__name__
        s += '%s(' % name
        s = _debug_args(expr.args, s, indent + len(name) + 1)

    elif isinstance(expr, AppliedUndef):
        name = expr.func.__name__
        s += name
        s = _debug_args(expr.args, s, indent + len(name) + 1)

    return s


def symbol_delete(sym):
    """Delete symbol.  This is useful if a symbol needs to be redefined
    with different assumptions."""

    state.symbols.pop(sym)


def symbol_map(name):

    return state.symbols.lookup(name)


# The following domain symbols are all SymPy symbols.
ssym = domainsymbol('s', complex=True)
tsym = domainsymbol('t', real=True)
f0sym = miscsymbol('f_0', real=True)
fsym = domainsymbol('f', real=True)
omegasym = domainsymbol('omega', real=True)
omega0sym = miscsymbol('omega_0', real=True)
tausym = miscsymbol('tau', real=True)
nusym = miscsymbol('nu', real=True)
Omegasym = domainsymbol('Omega', real=True)
Fsym = domainsymbol('F', real=True)
nsym = domainsymbol('n', integer=True)
ksym = domainsymbol('k', integer=True)
zsym = domainsymbol('z', real=False)

dt = miscsymbol('Delta_t', real=True, positive=True)
df = miscsymbol('Delta_f', real=True, positive=True)
eps = miscsymbol('epsilon42', real=True, positive=True)

pi = sym.pi
j = sym.I
oo = sym.oo
inf = sym.oo
one = sym.S.One

# This is required for expr('I') to work as j but it prevents having
# a current called I.
state.symbols.register('I', sym.I, 'misc')

try:
    from sympy.core.function import AppliedUndef
except:
    from sympy.function import AppliedUndef
