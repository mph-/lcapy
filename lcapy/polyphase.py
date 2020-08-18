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


class PolyphaseVoltageVector(PolyphaseVector):
    pass


class PolyphaseCurrentVector(PolyphaseVector):
    pass


class PolyphasePhase(object):

    def sequence(self):
        A = polyphase_decompose_matrix(len(self))
        return A * self    


class PolyphaseSequence(object):

    def phase(self):
        A = polyphase_compose_matrix(len(self))
        return A * self

    
class PhaseVoltageVector(PolyphaseVoltageVector, PolyphasePhase):

    def sequence(self):
        """Convert to sequence voltage vector."""
        return SequenceVoltageVector(super(PhaseVoltageVector, self).sequence())

    def line(self):
        """Convert to line voltage vector."""        
        D = phase_to_line_matrix(len(self))
        return LineVoltageVector(D * self)

    
class PhaseCurrentVector(PolyphaseCurrentVector, PolyphasePhase):

    def sequence(self):
        """Convert to sequence current vector."""        
        return SequenceCurrentVector(super(PhaseCurrentVector, self).sequence())    

    def line(self):
        """Convert to line current vector."""
        D = phase_to_line_matrix(len(self))
        return LineCurrentVector(D * self)
    

class SequenceVoltageVector(PolyphaseVoltageVector, PolyphaseSequence):

    def phase(self):
        """Convert to phase voltage vector."""                
        return PhaseVoltageVector(super(SequenceVoltageVector, self).phase())    

    def line(self):
        """Convert to line voltage vector."""        
        return self.phase().line()

    
class SequenceCurrentVector(PolyphaseCurrentVector, PolyphaseSequence):

    def phase(self):
        """Convert to phase current vector."""                        
        return PhaseCurrentVector(super(SequenceCurrentVector, self).phase())

    def line(self):
        """Convert to line current vector."""        
        return self.phase().line()
    
    
class LineVoltageVector(PolyphaseVoltageVector):
    """These are also known as phase to phase voltages."""
    pass


class LineCurrentVector(PolyphaseCurrentVector):
    """These are also known as phase to phase currents."""
    pass


def phase_to_line_matrix(N=3):

    a = Matrix.zeros(N)
    
    for row in range(N):
        a[row, row] = 1
        col = (row + 1) % N
        a[row, col] = -1
    return a


def polyphase_decompose_matrix(N=3, alpha=None):
    """Matrix to decompose vector of phase components into the
    sequence components.  The matrix dimension is `N` x `N`.  This
    matrix is equivalent to IDFTmatrix.  The transformation is only
    valid for positive frequency components; the complex conjugate is
    required for negative frequency components."""

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
    """Matrix to compose sequence components into a vector of phase
    components.  The matrix dimension is `N` x `N`.  This
    matrix is equivalent to DFTmatrix.  The transformation is only
    valid for positive frequency components; the complex conjugate is
    required for negative frequency components."""

    alpha = exp(-j * 2 * pi / N)

    a = Matrix.zeros(N)
    
    for row in range(N):
        for col in range(N):
            a[row, col] = alpha ** (row * col)
    return a


def polyphase_alpha(N):

    return exp(-j * 2 * pi / N)
