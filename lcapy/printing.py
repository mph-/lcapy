"""This module provides printing support.

Copyright 2014--2020 Michael Hayes, UCECE

"""

import re
from .config import latex_expr_map, pretty_expr_map, str_expr_map
from .config import functions, words, subscripts
from .latex import latex_str
from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.pretty.pretty import PrettyPrinter
import sympy as sym

__all__ = ('pretty', 'pprint', 'latex', 'print_str')

# Note, there are some magic hooks that are used for printing custom types:
#
# jupyter looks for methods called _repr_latex_ and _repr_pretty_.
#
# IPython looks for _repr_pretty
#
# sympy.latex() looks for methods called _latex

# LaTeX markup is nicer but it requires mathjax.

# FIXME: should import from parser
cpt_names = ('C', 'E', 'F', 'G', 'H', 'I', 'L', 'R', 'V', 'Y', 'Z', 'i', 'v')
cpt_name_pattern = re.compile(r"(%s)([\w']*)" % '|'.join(cpt_names))
sub_super_pattern = re.compile(r"([_\^]){([\w]+)}")
func_pattern = re.compile(r"\\operatorname{(.*)}")
word_name_pattern = re.compile(r"(%s)([\w']*)" % '|'.join(words))


def canonical_name(name):
    """Convert symbol name to canonical form for printing.

    R_{out} -> R_out
    R1 -> R_1
    XT2 -> XT_2
    Vbat -> V_bat
    alpha0 -> alpha_0
    """

    def foo(match):
        return match.group(1) + match.group(2)

    if not isinstance(name, str):
        return name

    # Convert R_{out} to R_out for SymPy to recognise.
    name = sub_super_pattern.sub(foo, name)

    # Convert v_Cfoo to v_C_foo, etc.   This is required for
    # state-space state variables.
    if (len(name) >= 3 and name[0:2] in ('v_', 'V_', 'i_', 'I_') and
        name[2] in ('L', 'C')):
        return name[0:2] + canonical_name(name[2:])
    
    if name.find('_') != -1:
        return name

    # Don't touch things like heaviside
    if name.lower() in words + functions:
        return name
    
    # Convert R1 to R_1, etc.
    match = cpt_name_pattern.match(name)
    if match:
        if match.groups()[1] == '':
            return name
        name = match.groups()[0] + '_' + match.groups()[1]
        return name

    if len(name) < 2:
        return name

    # Convert omega1 to omega_1, etc.
    match = word_name_pattern.match(name)
    if match:
        if match.groups()[1] == '':
            return name
        name = match.groups()[0] + '_' + match.groups()[1]
        return name
    
    # Convert i1 to i_1, etc.
    if name[1].isdigit:
        return name[0] + '_' + name[1:]

    # Convert irms to i_rms, etc.
    if name[1:].lower() in subscripts:
        return name[0] + '_' + name[1:]    

    return name


class LcapyStrPrinter(StrPrinter):

    def _print(self, expr, exp=None):

        from .expr import Expr
        if isinstance(expr, Expr):
            expr = expr.expr

        # Convert sym.I to j etc.
        try:            
            if expr in pretty_expr_map:        
                return str_expr_map[expr]
        except:
            pass
        return super(LcapyStrPrinter, self)._print(expr)

    def _print_Symbol(self, expr):
        # Do not canonicalise name since this is required for name
        # matching.  The only caller is the __str__ method for Expr.
        expr = sym.Symbol(expr.name)
        return super(LcapyStrPrinter, self)._print_Symbol(expr)    


class LcapyLatexPrinter(LatexPrinter):

    def _print(self, expr, exp=None):

        from .expr import Expr
        if isinstance(expr, Expr):        
            expr = expr.expr

        # Convert sym.I to j etc.
        try:
            if expr in latex_expr_map:
                return latex_expr_map[expr]
        except:
            pass

        if exp is None:
            return super(LcapyLatexPrinter, self)._print(expr)
        return super(LcapyLatexPrinter, self)._print(expr, exp=exp)        

    def _print_Piecewise(self, expr):

        if len(expr.args) > 1:
            return super(LcapyLatexPrinter, self)._print(expr)                    

        e, c = expr.args[0]
        return  r"%s \;\; \text{for}\: %s" % (self._print(e), self._print(c))

    def _print_Heaviside(self, expr, exp=None):

        tex = r"%s\left(%s\right)" % (latex_expr_map[sym.Heaviside],
                                      self._print(expr.args[0]))
        if exp:
            tex = r"\left(%s\right)^{%s}" % (tex, exp)
        return tex

    def _print_UnitImpulse(self, expr, exp=None):

        tex = r"\delta\left[%s\right]" % self._print(expr.args[0])
        if exp:
            tex = r"\delta\left[%s\right]^{%s}" % (tex, exp)
        return tex

    def _print_UnitStep(self, expr, exp=None):

        tex = r"u\left[%s\right]" % self._print(expr.args[0])
        if exp:
            tex = r"u\left[%s\right]^{%s}" % (tex, exp)
        return tex        

    def _print_symbol_name(self, name):

        name = canonical_name(name)
        
        parts = name.split('_')

        expr = sym.Symbol(parts[0])
        s = super(LcapyLatexPrinter, self)._print_Symbol(expr)

        if len(parts) == 1:
            return latex_str(s)

        if len(parts) == 2:
            return latex_str(s + '_{%s}' % parts[1])

        # Sympy cannot print a symbol name with a double subscript
        # using LaTeX.  This should be fixed in Sympy.
        # Need to convert v_C_1 to v_{C_{1}}        

        if len(parts) == 3:
            return latex_str(s + '_{%s_{%s}}' % (parts[1], parts[2]))

        raise ValueError('Cannot handle more than two subscripts for %s' % name)
    
    def _print_Symbol(self, expr):

        return self._print_symbol_name(expr.name)

    def _print_AppliedUndef(self, expr):

        s = self._print_symbol_name(expr.func.__name__)
        args = [str(self._print(arg)) for arg in expr.args]        
        return '%s(%s)' % (s, ','.join(args))

    
class LcapyPrettyPrinter(PrettyPrinter):

    def _print(self, expr, exp=None):

        from .expr import Expr
        if isinstance(expr, Expr):        
            expr = expr.expr

        try:            
            if expr in pretty_expr_map:
                return self._print_basestring(pretty_expr_map[expr])
        except:
            pass
        return super(LcapyPrettyPrinter, self)._print(expr)

    def _print_Symbol(self, expr):

        expr = sym.Symbol(canonical_name(expr.name))                
        parts = expr.name.split('_')        
        if len(parts) >= 2:
            # Due to unicode limitations, Sympy cannot print a symbol
            # name with a double subscript.  As a work-around combine
            # the subscripts.  Note, Sympy converts 'v_C1' into
            # 'v_C_1' so we need to clean up.
            expr.name = parts[0] + '_' + ''.join(parts[1:])
        s = super(LcapyPrettyPrinter, self)._print_Symbol(expr)
        return s

    def _print_Piecewise(self, expr):

        from sympy.printing.pretty.stringpict import prettyForm
        
        if len(expr.args) > 1:
            return super(LcapyPrettyPrinter, self)._print(expr)

        ec = expr.args[0]

        pform = self._print(ec.expr)
        pform = prettyForm(*pform.right('  for '))
        pform = prettyForm(*pform.right(self._print(ec.cond)))
        return pform

    def _print_Heaviside(self, expr):

        from sympy.printing.pretty.stringpict import prettyForm
        
        if self._use_unicode:
            pform = self._print(expr.args[0])
            pform = prettyForm(*pform.parens(left='(', right=')'))
            pform = prettyForm(*pform.left(pretty_expr_map[sym.Heaviside]))
            return pform
        else:
            return self._print_Function(expr)

    def _print_UnitImpulse(self, expr):

        from sympy.printing.pretty.stringpict import prettyForm
        from sympy.printing.pretty.pretty_symbology import greek_unicode
        
        if self._use_unicode:
            pform = self._print(expr.args[0])
            pform = prettyForm(*pform.parens(left='[', right=']'))
            pform = prettyForm(*pform.left(greek_unicode['delta']))
            return pform
        else:
            return self._print_Function(expr)

    def _print_UnitStep(self, expr):

        from sympy.printing.pretty.stringpict import prettyForm
        
        if self._use_unicode:
            pform = self._print(expr.args[0])
            pform = prettyForm(*pform.parens(left='[', right=']'))
            pform = prettyForm(*pform.left('u'))
            return pform
        else:
            return self._print_Function(expr)        
    
    
def print_str(expr):
    """Convert expression into a string."""
    
    return LcapyStrPrinter().doprint(expr)


def pretty(expr, **settings):
    """Pretty print an expression."""

    return LcapyPrettyPrinter(settings).doprint(expr)


def pprint(expr, **kwargs):
    """Pretty print an expression.

    If have non-interactive shell a latex string is returned."""

    import sys

    # If interactive use pretty, otherwise use latex
    if hasattr(sys, 'ps1'):
        print(pretty(expr, **kwargs))
    else:
        print(latex(expr, **kwargs))


def latex(expr, fold_frac_powers=False, fold_func_brackets=False,
    fold_short_frac=None, inv_trig_style="abbreviated",
    itex=False, ln_notation=False, long_frac_ratio=None,
    mat_delim="[", mat_str=None, mode="plain", mul_symbol=None,
    order=None, symbol_names=None):

    # This is mostly lifted from sympy/printing/latex.py when all we needed
    # was a hook...
    
    if symbol_names is None:
        symbol_names = {}

    settings = {
        'fold_frac_powers' : fold_frac_powers,
        'fold_func_brackets' : fold_func_brackets,
        'fold_short_frac' : fold_short_frac,
        'inv_trig_style' : inv_trig_style,
        'itex' : itex,
        'ln_notation' : ln_notation,
        'long_frac_ratio' : long_frac_ratio,
        'mat_delim' : mat_delim,
        'mat_str' : mat_str,
        'mode' : mode,
        'mul_symbol' : mul_symbol,
        'order' : order,
        'symbol_names' : symbol_names,
    }

    string = LcapyLatexPrinter(settings).doprint(expr)

    match = func_pattern.match(string)
    if match is not None:
        # v_1(t) -> \operatorname{v\_1}\left( t \right)
        # operatorname requires amsmath so switch to mathrm
        string = r'\mathrm{%s}' % match.groups()[0].replace('\\_', '_')

    return string


class PrintingConfig(object):

    def __init__(self):

        from sympy import init_printing
        init_printing(latex_printer=latex, pretty_printer=pretty, str_printer=print_str)

        self._abbreviate_units = False

    @property
    def abbreviate_units(self):

        return self._abbreviate_units

    @abbreviate_units.setter
    def abbreviate_units(self, val):    

        self._abbreviate_units = val
    
        # Print abbreviated units, V not volt
        sym.printing.str.StrPrinter._default_settings['abbrev'] = val
        

# See sympy/interactive/printing.py and IPython/core/formatters.py
# Also see hack at end of expr.py to support latex for Lcapy container
# types.




