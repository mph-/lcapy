import re
from .config import print_expr_map
from .latex import latex_str
from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.pretty.pretty import PrettyPrinter
import sympy as sym


__all__ = ('pretty', 'pprint', 'latex', 'print_str')

# Note, jupyter looks for methods called _repr_latex_ and
# _repr_pretty_.  LaTeX markup is nicer but it requires mathjax.

cpt_names = ('C', 'E', 'F', 'G', 'H', 'I', 'L', 'R', 'V', 'Y', 'Z')
cpt_name_pattern = re.compile(r"(%s)([\w']*)" % '|'.join(cpt_names))
sub_super_pattern = re.compile(r"([_\^]){([\w]+)}")
func_pattern = re.compile(r"\\operatorname{(.*)}")


def canonical_name(name):
    """Convert symbol name to canonical form for printing.

    R_{out} -> R_out
    R1 -> R_1
    XT2 -> XT_2
    Vbat -> V_bat
    """

    def foo(match):
        return match.group(1) + match.group(2)

    if not isinstance(name, str):
        return name

    # Convert R_{out} to R_out for SymPy to recognise.
    name = sub_super_pattern.sub(foo, name)

    if name.find('_') != -1:
        return name

    # TODO, add other special function names
    if name in ('Heaviside', 'DiracDelta'):
        return name
    
    # Rewrite R1 as R_1, etc.
    match = cpt_name_pattern.match(name)
    if match:
        if match.groups()[1] == '':
            return name
        name = match.groups()[0] + '_' + match.groups()[1]
        return name

    return name


class LcapyStrPrinter(StrPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr
            if expr in print_expr_map:
                return print_expr_map[expr]        
        return super(LcapyStrPrinter, self)._print(expr)

    def _print_Symbol(self, expr):
        expr = sym.Symbol(canonical_name(expr.name))                
        return super(LcapyStrPrinter, self)._print_Symbol(expr)    


class LcapyLatexPrinter(LatexPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr            
            if expr in print_expr_map:
                return print_expr_map[expr]
        return super(LcapyLatexPrinter, self)._print(expr)

    def _print_Symbol(self, expr):
        expr = sym.Symbol(canonical_name(expr.name))                
        return super(LcapyLatexPrinter, self)._print_Symbol(expr)    


class LcapyPrettyPrinter(PrettyPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr
            if expr in print_expr_map:
                return self._print_basestring(print_expr_map[expr])
        return super(LcapyPrettyPrinter, self)._print(expr)

    def _print_Symbol(self, expr):
        expr = sym.Symbol(canonical_name(expr.name))                
        return super(LcapyPrettyPrinter, self)._print_Symbol(expr)


def print_str(expr):
    """Convert expression into a string."""
    
    return LcapyStrPrinter().doprint(expr)


def latex(expr, **settings):
    """Convert expression into a LaTeX string."""
    
    string = LcapyLatexPrinter(settings).doprint(expr)

    match = func_pattern.match(string)
    if match is not None:
        # v_1(t) -> \operatorname{v\_1}\left( t \right)
        # operatorname requires amsmath so switch to mathrm
        string = r'\mathrm{%s}' % match.groups()[0].replace('\\_', '_')

    return latex_str(string)


def pretty(expr, **settings):
    """Pretty print an expression."""
    
    return LcapyPrettyPrinter(settings).doprint(expr)


def pprint(expr):
    """Pretty print an expression.

    If have non-interactive shell a latex string is returned."""

    import sys

    # If interactive use pretty, otherwise use latex
    if hasattr(sys, 'ps1'):
        print(pretty(expr))
    else:
        print(latex(expr))

