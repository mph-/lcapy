"""
This module implements the Vector class (a wrapper for Matrix).

Copyright 2020--2021 Michael Hayes, UCECE
"""

from __future__ import division
from .matrix import Matrix
from .sym import sympify
from .expr import expr


class Vector(Matrix):

    def __new__(cls, *args, **assumptions):

        if len(args) == 2:
            return super(Vector, cls).__new__(cls, (expr(args[0], **assumptions).expr, expr(args[1], **assumptions).expr))

        args = [expr(arg, **assumptions).sympy for arg in args[0]]

        return super(Vector, cls).__new__(cls, args)
