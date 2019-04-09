import sympy as sym
from copy import copy
from .printing import pprint, latex
from .expr import Expr
from .laplace import laplace_transform, inverse_laplace_transform

def msympify(expr):
    # If do nothing, will get a problem with matrices that
    # have mixed data types, e.g., A matrix.

    if isinstance(expr, Expr):
        # Bye bye lcapy type information...
        return expr.expr
    return expr
        
        
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
    # prevent this by defining _simpify to do nothing instead of
    # sym.sympify.
    
    # _typewrap = sExpr

    _sympify = staticmethod(msympify)
    
    def __getitem__(self, key):

        item = super(Matrix, self).__getitem__(key)

        # The following line is to handle slicing used
        # by latex method.
        if isinstance(item, sym.Matrix):
            return item

        if hasattr(self, '_typewrap'):
            return self._typewrap(item)

        return item

    def pprint(self):

        return pprint(self)

    def latex(self):

        return latex(self)

    def _reformat(self, methodname):
        """Helper method for reformatting expression."""

        new = copy(self)

        for i in range(self.rows):
            for j in range(self.cols):
                method = getattr(self[i, j], methodname)
                new[i, j] = method()

        return new

    def laplace(self):

        def lt(expr):
            from .sym import ssym, tsym
            return laplace_transform(expr, tsym, ssym)
        
        return self.applyfunc(lt)

    def inverse_laplace(self, **assumptions):

        def lt(expr):
            from .sym import ssym, tsym
            return inverse_laplace_transform(expr, ssym, tsym, **assumptions)
        
        return self.applyfunc(lt)    

    
    def canonical(self):

        return self._reformat('canonical')

    def general(self):

        return self._reformat('general')

    def mixedfrac(self):

        return self._reformat('mixedfrac')

    def partfrac(self):

        return self._reformat('partfrac')

    def timeconst(self):

        return self._reformat('timeconst')    

    def ZPK(self):

        return self._reformat('ZPK')

    # TODO. There is probably a cunning way to automatically handle
    # the following.

    def inv(self):

        return self.__class__(sym.Matrix(self).inv())

    def det(self):

        return sym.Matrix(self).det()

    def simplify(self):
        # TODO
        return self

#from lcapy.sexpr import sExpr

