import sympy as sym
from .expr import expr


class SystemEquations(object):

    def __init__(self, A, b, y):
        """Represent a system of equations as A y = b."""

        self.A = A
        self.b = sym.Matrix(b)
        self.y = sym.Matrix(y)

    @property
    def Ainv(self):

        if hasattr(self, '_Ainv'):
            return self._Ainv
        self._Ainv = self.A.inv()
        return self._Ainv

    def format(self, form='A y = b', invert=False):
        """Forms can be:
         A y = b
         b = A y
         Ainv b = y
         y = Ainv b

        If `invert` is True, the A matrix is inverted."""

        if form == 'A y = b':
            return expr(sym.Eq(sym.MatMul(self.A, self.y), self.b), evaluate=False)

        elif form == 'b = A y':
            return expr(sym.Eq(self.b, sym.MatMul(self.A, self.y)), evaluate=False)

        elif form in ('y = Ainv b', 'default'):
            if invert:
                return expr(sym.Eq(self.y, sym.MatMul(self.Ainv, self.b),
                                   evaluate=False))

            return expr(sym.Eq(self.y, sym.MatMul(sym.Pow(self.A, -1), self.b),
                               evaluate=False))

        elif form == 'Ainv b = y':
            if invert:
                return expr(sym.Eq(sym.MatMul(self.Ainv, self.b), self.y,
                                   evaluate=False))

            return expr(sym.Eq(sym.MatMul(sym.Pow(self.A, -1), self.b), self.y,
                               evaluate=False))
        else:
            raise ValueError('Unknown form %s' % form)
