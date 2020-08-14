"""
This module implements a wrapper for the SymPy Matrix class.

Copyright 2019--2020 Michael Hayes, UCECE
"""

import sympy as sym
from copy import copy
from .sym import simplify
from .printing import pprint, latex, pretty


def msympify(expr):
    # If do nothing, will get a problem with matrices that
    # have mixed data types, e.g., A matrix.

    if isinstance(expr, Expr):
        # Bye bye Lcapy type information...
        return expr.expr
    return sym.sympify(expr)
        
        
class Matrix(sym.Matrix):

    # Unlike numpy.ndarray, the sympy.Matrix runs all the elements
    # through _sympify, creating sympy objects and thus losing the
    # original type information and associated methods.  As a hack, we
    # try to wrap elements when they are read using __getitem__.  This
    # assumes that all the elements have the same type.  This is not
    # the case for A, B, G, and H matrices.  This could he handled by
    # having another matrix to specify the type for each element.

    # What's worse, is that calling _sympify on each element creates
    # different variables than what we are expecting.  For example,
    # the sExpr s looks the same but gets different attributes.  We
    # prevent this by defining _sympify.
    
    # _typewrap = sExpr

    _sympify = staticmethod(msympify)
    
    def __getitem__(self, key):

        item = super(Matrix, self).__getitem__(key)

        # The following line is to handle slicing used
        # by the latex method.
        if isinstance(item, sym.Matrix):
            return item

        if hasattr(self, '_typewrap'):
            return self._typewrap(item)

        return item

    def _repr_pretty_(self, p, cycle):
        """This is used by jupyter notebooks to display an expression using
        unicode."""

        p.text(pretty(self))

    def _repr_latex_(self):
        """This is used by jupyter notebooks to display an expression using
        LaTeX markup.  However, this requires mathjax.  If this method
        is not defined, jupyter falls back on _repr__pretty_ which
        outputs unicode."""
        import pdb; pdb.set_trace()
        return '$$' + latex(self) + '$$'

    def pprint(self):

        return pprint(self)

    def latex(self, **kwargs):

        return latex(self, **kwargs)

    def canonical(self):
        
        return self

    # TODO. There is probably a cunning way to automatically handle
    # the following.

    def inv(self):

        return self.__class__(sym.Matrix(self).inv())

    def det(self):

        return sym.Matrix(self).det()

    def simplify(self):
        
        return self.applyfunc(simplify)

    @property
    def symbols(self):

        symbols = {}
        for elt in self:
            symbols.update(expr(elt).symbols)
        return symbols        
    

def matrix(mat):
    """Create Lcapy Matrix from a SymPy Matrix.

    If a t symbol is found in an element a tMatrix object is created.
    If a s symbol is found in an element an sMatrix object is created.

    """

    from .sym import tsym, ssym
    from .smatrix import sMatrix
    from .tmatrix import tMatrix
    
    elt = mat[0]
    try:
        elt = elt[0]
    except:
        pass    
    
    if elt.has(tsym):
        return tMatrix(mat)
    elif elt.has(ssym):
        return sMatrix(mat)        
    else:
        return mat

from .expr import Expr, expr
