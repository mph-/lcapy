"""
This module implements a quantity for the SymPy Matrix class for vectors.

Copyright 2020 Michael Hayes, UCECE
"""


from __future__ import division
from .matrix import Matrix
from .sym import sympify
from .expr import expr

class Vector(Matrix):

    def __new__(cls, *args):

        if len(args) == 2:
            return super(Vector, cls).__new__(cls, (expr(args[0]).expr, expr(args[1]).expr))

        args = [expr(arg).expr for arg in args[0]]

        return super(Vector, cls).__new__(cls, args)


