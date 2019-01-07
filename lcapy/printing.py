import sympy as sym
import re

from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.pretty.pretty import PrettyPrinter

func_pattern = re.compile(r"\\operatorname{(.*)}")


class LcapyStrPrinter(StrPrinter):

    def _print(self, expr):

        if hasattr(expr, 'expr'):
            expr = expr.expr

        if expr == sym.I:
            return "j"
        return super(LcapyStrPrinter, self)._print(expr)


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


def pprint(expr):

    # If interactive use pretty, otherwise use latex
    if hasattr(sys, 'ps1'):
        print(pretty(expr))
    else:
        print(latex_str(latex(expr)))

        
def print_str(expr):

    return LcapyStrPrinter().doprint(expr)
    
    
