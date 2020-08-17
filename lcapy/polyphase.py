"""This module provides preliminary support for polyphase systems. 

Copyright 2020 Michael Hayes, UCECE

"""

from .matrix import Matrix
from .sym import j, pi
from .functions import exp
from .expr import expr
from .phasor import Vphasor, Iphasor
from .vector import Vector


class PolyphaseVector(Vector):
    pass


class PolyphaseVoltage(PolyphaseVector):
    pass


class PolyphaseCurrent(PolyphaseVector):
    pass


class PolyphaseUnbalanced(object):

    def sequence(self):

        A = polyphase_decompose_matrix(len(self))
        return A * self    


class PolyphaseSequence(object):

    def unbalanced(self):

        A = polyphase_compose_matrix(len(self))
        return A * self

    
class PolyphaseUnbalancedVoltage(PolyphaseVoltage, PolyphaseUnbalanced):

    def sequence(self):
        return PolyphaseSequenceVoltage(super(PolyphaseUnbalancedVoltage, self).sequence())


class PolyphaseUnbalancedCurrent(PolyphaseCurrent, PolyphaseUnbalanced):

    def sequence(self):
        return PolyphaseSequenceCurrent(super(PolyphaseUnbalancedCurrent, self).sequence())    


class PolyphaseSequenceVoltage(PolyphaseVoltage, PolyphaseSequence):

    def unbalanced(self):
        return PolyphaseUnbalancedVoltage(super(PolyphaseSequenceVoltage, self).unbalanced())    


class PolyphaseSequenceCurrent(PolyphaseCurrent, PolyphaseSequence):

    def unbalanced(self):
        return PolyphaseUnbalancedCurrent(super(PolyphaseSequenceCurrent, self).unbalanced())        



def polyphase_decompose_matrix(N=3, alpha=None):
    """Matrix to decompose vector of unbalanced polyphase phasors into the
    sequence components.  The matrix dimension is `N` x `N`.  This
    matrix is equivalent to IDFTmatrix.

    """

    if alpha is None:
        alpha = exp(j * 2 * pi / N)
    else:
        alpha = expr(alpha)

    a = Matrix.zeros(N)
    
    for row in range(N):
        for col in range(N):
            a[row, col] = alpha ** (row * col)
    return a / N


def polyphase_compose_matrix(N=3):
    """Matrix to compose sequence components into a vector of unbalanced
    polyphase phasors.  The matrix dimension is `N` x `N`.  This
    matrix is equivalent to DFTmatrix."""

    alpha = exp(-j * 2 * pi / N)

    a = Matrix.zeros(N)
    
    for row in range(N):
        for col in range(N):
            a[row, col] = alpha ** (row * col)
    return a
