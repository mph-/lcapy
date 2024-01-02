"""This module provides the SystemEquations class.

Copyright 2023--2024 Michael Hayes, UCECE

"""


import sympy as sym
from .expr import equation


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
            return equation(sym.MatMul(self.A, self.y), self.b)

        elif form == 'b = A y':
            return equation(self.b, sym.MatMul(self.A, self.y))

        elif form in ('y = Ainv b', 'default'):
            if invert:
                return equation(self.y, sym.MatMul(self.Ainv, self.b))

            # Shown as A^{-1}; inverse not computed
            Ainv = sym.Pow(self.A, -1)
            return equation(self.y, sym.MatMul(Ainv, self.b))

        elif form == 'Ainv b = y':
            if invert:
                return equation(sym.MatMul(self.Ainv, self.b), self.y,
                                evaluate=False)

            # Shown as A^{-1}; inverse not computed
            Ainv = sym.Pow(self.A, -1)
            return equation(sym.MatMul(Ainv, self.b), self.y,
                            evaluate=False)
        else:
            raise ValueError('Unknown form %s' % form)
