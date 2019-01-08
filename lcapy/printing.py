import sympy as sym
import re
from .config import print_expr_map
from .latex import latex_str
from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.pretty.pretty import PrettyPrinter

__all__ = ('pretty', 'pprint', 'latex', 'print_str')


func_pattern = re.compile(r"\\operatorname{(.*)}")


class LcapyStrPrinter(StrPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr
            if expr in print_expr_map:
                return print_expr_map[expr]        
        return super(LcapyStrPrinter, self)._print(expr)


class LcapyLatexPrinter(LatexPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr
            if expr in print_expr_map:
                return print_expr_map[expr]
        return super(LcapyLatexPrinter, self)._print(expr)


class LcapyPrettyPrinter(PrettyPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr
            if expr in print_expr_map:
                return self._print_basestring(print_expr_map[expr])
        return super(LcapyPrettyPrinter, self)._print(expr)


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
