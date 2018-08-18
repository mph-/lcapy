from sympy.parsing.sympy_parser import parse_expr, auto_number, rationalize
try:
    from sympy.parsing.sympy_parser import NUMBER, STRING, NAME, OP        
except:
    from sympy.parsing.sympy_tokenize import NUMBER, STRING, NAME, OP
    
from sympy import Basic, Symbol, Expr
from sympy.core.function import AppliedUndef
import sympy as sym
import re

global_dict = {}
exec('from sympy import *', global_dict)
global_ignore = ('C', 'O', 'S', 'N', 'E', 'E1', 'Q')
for symbol in global_ignore:
    global_dict.pop(symbol)
# delta gets printed as DiracDelta; could override
global_dict['delta'] = global_dict['DiracDelta']
global_dict['step'] = global_dict['Heaviside']
global_dict['u'] = global_dict['Heaviside']

cpt_names = ('C', 'E', 'F', 'G', 'H', 'I', 'L', 'R', 'V', 'Y', 'Z')
cpt_name_pattern = re.compile(r"(%s)([\w']*)" % '|'.join(cpt_names))

sub_super_pattern = re.compile(r"([_\^]){([\w]+)}")

def canonical_name(name):

    def foo(match):
        return match.group(1) + match.group(2)

    if not isinstance(name, str):
        return name

    # Convert R_{out} to R_out for sympy to recognise.
    name = sub_super_pattern.sub(foo, name)

    if name.find('_') != -1:
        return name

    # Rewrite R1 as R_1, etc.
    match = cpt_name_pattern.match(name)
    if match:
        if match.groups()[1] == '':
            return name
        name = match.groups()[0] + '_' + match.groups()[1]
        return name

    return name

def symbols_find(arg):
    """Return list of symbols in arg.  No symbols are cached."""

    symbols = []

    def find_symbol(tokens, local_dict, global_dict):
        
        for tok in tokens:
            tokNum, tokVal = tok
            if tokNum == NAME:
                name = tokVal
                if name == 'j':
                    name = 'I'
                if name not in local_dict and name not in global_dict:
                    symbols.append(name)
        return ([(NUMBER, '0')])

    if isinstance(arg, str):
        parse_expr(arg, transformations=(find_symbol, ), 
                   global_dict=global_dict, local_dict={}, evaluate=False)
        
        return symbols

    # Hack
    if hasattr(arg, 'expr'):
        arg = arg.expr

    if not isinstance(arg, (Symbol, Expr, AppliedUndef)):
        return []
    return [repr(symbol) for symbol in arg.atoms(Symbol, AppliedUndef)]

def parse(string, symbols={}, evaluate=True, local_dict={}, **assumptions):
    """Handle arbitrary strings that may refer to multiple symbols."""

    cache = assumptions.pop('cache', True)

    def auto_symbol(tokens, local_dict, global_dict):
        """Inserts calls to ``Symbol`` or ``Function`` for undefined variables/functions."""
        result = []
        prevTok = (None, None)

        tokens.append((None, None))  # so zip traverses all tokens
        for tok, nextTok in zip(tokens, tokens[1:]):
            tokNum, tokVal = tok
            nextTokNum, nextTokVal = nextTok
            if tokNum == NAME:
                name = tokVal
                if name == 'j':
                    name = 'I'

                if name in global_dict:

                    obj = global_dict[name]
                    if isinstance(obj, (Basic, type)):
                        result.append((NAME, name))
                        continue

                    if callable(obj):
                        result.append((NAME, name))
                        continue

                name = canonical_name(str(name))

                if name in local_dict:
                    # print('Found %s' % name)
                    # Could check assumptions.
                    result.append((NAME, name))
                    continue

                # Automatically add Function.  We ignore the assumptions.
                # These could be supported by modifying fourier.py/laplace.py
                # to propagate assumptions when converting V(s) to v(t), etc.
                if nextTokVal == '(':
                    result.extend([(NAME, 'Function'),
                                   (OP, '('), (NAME, repr(name)), (OP, ')')])
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

            prevTok = (tokNum, tokVal)

        return result


    s = parse_expr(string, transformations=(auto_symbol, auto_number,
                                            rationalize), 
                   global_dict=global_dict, local_dict=local_dict,
                   evaluate=evaluate)
    if not cache:
        return s

    # Look for newly defined symbols/functions.
    for symbol in s.atoms(Symbol, AppliedUndef):
        name = repr(symbol)
        if name not in symbols:
            symbols[name] = symbol

    return s


def sympify1(arg, symbols={}, evaluate=True, **assumptions):
    """Create a sympy expression.

    The purpose of this function is to head sympy off at the pass and
    apply the defined assumptions.

    """

    if hasattr(arg, 'expr'):
        return arg.expr

    if isinstance(arg, (Symbol, Expr)):
        return arg

    # Why doesn't sympy do this?
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
        # Handle arbitrary strings that may refer to multiple symbols.
        return parse(arg, symbols, evaluate=evaluate,
                     local_dict=symbols, **assumptions)

    return sym.sympify(arg, rational=True, locals=symbols, 
                       evaluate=evaluate)

def test():
    symbols = {}
    s1 = sympify1('5 * E1 + a', symbols)
    s2 = sympify1('5 * R1 + a', symbols, real=True)
    print(symbols['R_1'].assumptions0)
    s3 = sympify1('5 * R1 + a', symbols, positive=True)
    print(symbols1['R_1'].assumptions0)

