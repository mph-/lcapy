"""This module defines the StateSpace class for representing a linear
time-invariant system as a state-space model.

Copyright 2019-2021 Michael Hayes, UCECE

"""

from .smatrix import LaplaceDomainMatrix
from .tmatrix import TimeDomainMatrix
from .sym import ssym
from .texpr import t
import sympy as sym


class StateSpace(object):

    def __init__(self, A, B, C, D, u, y, x, x0):

        self._A = A
        self._B = B
        self._C = C
        self._D = D
        self._u = u
        self._x = x
        self._x0 = x0        
        self._y = y

        self.dotx = TimeDomainMatrix([sym.Derivative(x1, t) for x1 in x])
    
    def state_equations(self):
        """System of first-order differential state equations:

        dotx = A x + B u

        where x is the state vector and u is the input vector.
        """
        
        return expr(sym.Eq(self.dotx, sym.MatAdd(sym.MatMul(self._A, self._x),
                                                 sym.MatMul(self._B, self._u)),
                           evaluate=False))

    def output_equations(self):
        """System of output equations:

        y = C x + D u

        where y is the output vector, x is the state vector and u is
        the input vector.

        """
        
        return expr(sym.Eq(self._y, sym.MatAdd(sym.MatMul(self._C, self._x),
                                               sym.MatMul(self._D, self._u)),
                           evaluate=False))

    @property
    def x(self):
        """State variable vector."""
        return self._x

    @property
    def x0(self):
        """State variable initial value vector."""
        return self._x0

    @property
    def u(self):
        """Input vector."""
        return self._u
    
    @property
    def y(self):
        """Output vector."""
        return self._y    

    @property
    def A(self):
        """State matrix."""
        return self._A

    @property
    def B(self):
        """Input matrix."""
        return self._B

    @property
    def C(self):
        """Output matrix."""
        return self._C

    @property
    def D(self):
        """Feed-through matrix."""
        return self._D        

    @property
    def Phi(self):
        """s-domain state transition matrix."""

        M = LaplaceDomainMatrix(sym.eye(len(self._x)) * ssym - self._A)
        return LaplaceDomainMatrix(M.inv().canonical())

    @property
    def phi(self):
        """State transition matrix."""        
        return TimeDomainMatrix(self.Phi.inverse_laplace(causal=True))
        
    @property
    def U(self):
        """Laplace transform of input vector."""
        return LaplaceDomainMatrix(self._u.laplace())

    @property
    def X(self):
        """Laplace transform of state-variable vector."""        
        return LaplaceDomainMatrix(self._x.laplace())

    @property
    def Y(self):
        """Laplace transform of output vector."""        
        return LaplaceDomainMatrix(self._y.laplace())

    @property
    def H(self):
        """X(s) / U(s)"""

        return LaplaceDomainMatrix(self.Phi * self._B).canonical()

    @property
    def h(self):
        return TimeDomainMatrix(self.H.inverse_laplace(causal=True))

    @property
    def G(self):
        """System transfer functions."""

        return LaplaceDomainMatrix(self._C * self.H + self._D).canonical()

    @property
    def g(self):
        """System impulse responses."""        
        return TimeDomainMatrix(self.G.inverse_laplace(causal=True))
    
    def characteristic_polynomial(self):
        """Characteristic polynomial (aka system polynomial).

        `lambda(s) = |s * I - A|`"""

        M = Matrix(sym.eye(len(self._x)) * ssym - self._A)        
        return LaplaceDomainExpression(M.det()).simplify()

    @property
    def P(self):
        """Characteristic polynomial (aka system polynomial).

        `lambda(s) = |s * I - A|`"""        

        return self.characteristic_polynomial().canonical()

    @property        
    def eigenvalues_dict(self):
        """Dictionary of eigenvalues, the roots of the characteristic
        polynomial (equivalent to the poles of Phi(s)).  The
        dictionary values are the multiplicity of the eigenvalues.

        For a list of eigenvalues use eigenvalues."""        

        return self.characteristic_polynomial().roots()
        
    @property        
    def eigenvalues(self):
        """List of eigenvalues, the roots of the characteristic polynomial
        (equivalent to the poles of Phi(s))."""
        
        roots = self.eigenvalues_dict
        e = []

        # Replicate duplicated eigenvalues and return as a list.
        for v, n in roots.items():
            for m in range(n):
                e.append(v)
        return ExprList(e)

    @property    
    def Lambda(self):
        """Diagonal matrix of eigenvalues."""

        # Perhaps faster to use diagonalize
        # E, L = self.A.diagonalize()
        # return L
        
        e = self.eigenvalues
        return LaplaceDomainMatrix(sym.diag(*e))

    @property        
    def eigenvectors(self):
        """List of tuples (eigenvalue, multiplicity of eigenvalue,
        basis of the eigenspace) of A."""
        
        return self._A.eigenvects()
    
    @property    
    def M(self):
        """Modal matrix (eigenvectors of A)."""

        E, L = self._A.diagonalize()
        
        return LaplaceDomainMatrix(E)
    
    
from .symbols import t, s
from .expr import ExprList
from .sexpr import LaplaceDomainExpression
