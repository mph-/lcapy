from __future__ import division
from .matrix import Matrix
from .sym import sympify


class Vector(Matrix):

    def __new__(cls, *args):

        args = [sympify(arg) for arg in args]

        if len(args) == 2:
            return super(Vector, cls).__new__(cls, (args[0], args[1]))

        return super(Vector, cls).__new__(cls, *args)


