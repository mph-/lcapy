"""This module contains functions for handling symbols and in particular
for integrating with SymPy.

Copyright 2014--2020 Michael Hayes, UCECE

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
from .simplify import simplify_dirac_delta, simplify_heaviside

__all__ = ('symsymbol', 'sympify', 'simplify', 'symbol_delete')


global_dict = {}
exec('from sympy import *', global_dict)

for _alias, _name in aliases.items():
    global_dict[_alias] = global_dict[_name]

global_dict['abs'] = sym.Abs

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

    from .expr import Expr as LExpr
    if isinstance(arg, LExpr):
        arg = arg.expr
    
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
    for symbol in s.atoms(Symbol):
        name = symbol_name(symbol)
        if name not in symbols:
            if False:
                print("Adding symbol '%s'" % name)
                
            symbols[name] = symbol

    return s


def sympify1(arg, symbols=None, evaluate=True, override=False, **assumptions):
    """Create a SymPy expression.

    The purpose of this function is to head SymPy off at the pass and
    apply the defined assumptions.

    If `evaluate` is True, the expression is evaluated.

    If `override` is True, then create new symbol(s) even if
    previously defined by SymPy.
    """

    if symbols is None:
        symbols = {}

    if isinstance(arg, (Symbol, Expr)):
        return arg

    # Why doesn't SymPy do this?
    if isinstance(arg, complex):
        re = sym.sympify(str(arg.real), rational=True, evaluate=evaluate)
        im = sym.sympify(str(arg.imag), rational=True, evaluate=evaluate)
        if im == 1.0:
            arg = re + sym.I
        else:
            arg = re + sym.I * im
        return arg

    if isinstance(arg, float):
        # Note, need to convert to string to achieve a rational
        # representation.
        return sym.sympify(str(arg), rational=True, evaluate=evaluate)
        
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

    return sym.sympify(arg, rational=True, locals=symbols, 
                       evaluate=evaluate)


def sympify(expr, evaluate=True, override=False, **assumptions):
    """Create a SymPy expression.

    By default, symbols are assumed to be positive if no assumptions
    are explicitly defined.

    Note, this will not modify previously defined symbols with the
    same name.  Thus you cannot change the assumptions.

    If `evaluate` is True, the expression is evaluated.

    If `override` is True, then create new symbol(s) even if
    previously defined by SymPy.
    """
    
    if assumptions == {}:
        assumptions['positive'] = True
        # Note this implies that imag is False.   Also note that all
        # reals are considered complex (but with a zero imag part).

    elif 'positive' in assumptions:
        if not assumptions['positive']:
            assumptions.pop('positive')
        
    return sympify1(expr, state.context.symbols, evaluate, override,
                    **assumptions)


def symsymbol(name, **assumptions):
    """Create a SymPy symbol.

    By default, symbols are assumed to be positive unless real is
    defined.

    """
    return sympify(name, override=True, **assumptions)


def symsimplify(expr):
    """Simplify a SymPy expression.  This is a hack to work around
    problems with SymPy's simplify API."""

    # Handle Matrix types
    if hasattr(expr, 'applyfunc'):
        return expr.applyfunc(lambda x: symsimplify(x))

    if expr.has(sym.DiracDelta):
        expr = simplify_dirac_delta(expr)
    if expr.has(sym.Heaviside):
        expr = simplify_heaviside(expr)

    # This gets expanded into piecewise...
    if expr.has(sym.sign):
        # Could replace sign with dummy function, simplify, and then restore...
        # Should replace 1 + sign(t) with 2 * Heaviside(t)
        return expr
    
    try:
        if expr.is_Function and expr.func in (sym.Heaviside, sym.DiracDelta):
            return expr
    except:
        pass

    expr = sym.simplify(expr)

    return expr


def simplify(expr):
    """Simplify an Lcapy expression.  This is not straightforward, see
    sympy.simplify."""

    try:
        return expr.simplify()
    except:
        pass

    from .expr import Expr as LExpr
    if isinstance(expr, LExpr):
        expr = expr.expr

    return symsimplify(expr)


def is_sympy(expr):
    return isinstance(expr, (Symbol, Expr, AppliedUndef))


def symdebug(expr, s='', indent=0):

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
    
    state.context.symbols.pop(sym)
    

def symbol_map(name):

    new = name
    if not isinstance(name, str):
        name = str(name)
    
    # Replace symbol names with symbol definitions to
    # avoid problems with real or positive attributes.
    if name in state.context.symbols:
        new = state.context.symbols[name]
    elif name in state.global_context.symbols:
        new = state.global_context.symbols[name]            
    else:
        # Perhaps have symbol defined using sympy?
        pass
    return new


# The following are all SymPy symbols.
ssym = symsymbol('s', complex=True)
tsym = symsymbol('t', real=True)
fsym = symsymbol('f', real=True)
omegasym = symsymbol('omega', real=True)
omega0sym = symsymbol('omega_0', real=True)
tausym = symsymbol('tau', real=True)
nusym = symsymbol('nu', real=True)

pi = sym.pi
j = sym.I
oo = sym.oo
inf = sym.oo
one = sym.S.One

try:
    from sympy.core.function import AppliedUndef
except:
    from sympy.function import AppliedUndef    
