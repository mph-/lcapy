"""This module defines the StateSpace class for representing a linear
time-invariant system as a state-space model.

Copyright 2019-2021 Michael Hayes, UCECE

"""

from .matrix import Matrix
from .smatrix import LaplaceDomainMatrix
from .tmatrix import TimeDomainMatrix
from .sym import ssym
from .texpr import t, texpr
from .expr import expr
import sympy as sym


class StateSpace(object):

    def __init__(self, A, B, C, D, u=None, y=None, x=None, x0=None):
        """A is Nx x Nx state matrix
        B is Nx x Nu input matrix
        C is Ny x Nx output matrix
        D is Ny x Nu feedthrough matrix

        u is Nu x 1 input vector
        x is Nx x 1 state vector
        x0 is Nx x 1 state initial vector
        y is Ny x 1 state output vector
        """

        if not isinstance(A, Matrix):
            raise ValueError('A not matrix')
        if not isinstance(B, Matrix):
            raise ValueError('B not matrix')
        if not isinstance(C, Matrix):
            raise ValueError('C not matrix')
        if not isinstance(D, Matrix):
            raise ValueError('D not matrix')        

        # Number of state variables (the system order).
        Nx = A.shape[0]

        # Number of inputs
        Nu = B.shape[1]

        # Number of outputs
        Ny = C.shape[0]

        if A.shape[0] != A.shape[1]:
            raise ValueError('A matrix not square')

        if B.shape[0] != Nx:
            raise ValueError('B matrix has wrong dimension')
        if C.shape[1] != Nx:
            raise ValueError('C matrix has wrong dimension')        
        if (D.shape[0] != Ny) or (D.shape[1] != Nu):
            raise ValueError('D matrix has wrong dimension')                

        if u is None:
            u = TimeDomainMatrix([texpr('u_%d(t)' % n) for n in range(Nu)])

        if x is None:
            x = TimeDomainMatrix([texpr('x_%d(t)' % n) for n in range(Nx)])

        if x0 is None:
            x0 = x * 0

        if y is None:
            y = TimeDomainMatrix([texpr('y_%d(t)' % n) for n in range(Ny)])
            
        if x.shape[0] != Nx:
            raise ValueError('x vector has wrong dimension')
        if x0.shape[0] != Nx:
            raise ValueError('x0 vector has wrong dimension')
        if u.shape[0] != Nu:
            raise ValueError('u vector has wrong dimension')
        if y.shape[0] != Ny:
            raise ValueError('y vector has wrong dimension')                
            
        self._A = A
        self._B = B
        self._C = C
        self._D = D
        self._u = u
        self._x = x
        self._x0 = x0        
        self._y = y

        self.dotx = TimeDomainMatrix([sym.Derivative(x1, t) for x1 in x])

    @classmethod
    def from_transfer_function_coeffs(cls, b, a, method='CCF'):
        """Create state-space representation from transfer function
        specified with numerator and denominator coefficients.

         Note, state-space representations are not unique and are
        determined by the `method` argument.  Currently this can be
        'CCF' for the canonical controllable form.

        """        

        if method != 'CCF':
            raise ValueError('Only CCF method currently supported')

        b = list(b)
        a = list(a)
        Nb = len(b)
        Na = len(a)

        a0 = a[0]
        if a0 != 1:
            a = [ax / a0 for ax in a]
            b = [bx / a0 for bx in b]

        if Na > Nb:
            b = [0] * (Na - Nb) + b
        if Nb > Na:
            a = [0] * (Nb - Na) + a            

        Nx = Na - 1
        Nu = 1
        Ny = 1

        A = Matrix.zeros(Nx, Nx)
        B = Matrix.zeros(Nx, 1)
        C = Matrix.zeros(1, Nx)
        D = Matrix.zeros(1, 1)

        D[0, 0] = b[0]
        for n in range(Nx):
            C[0, n] = b[Nx - n] - a[Nx - n] * b[0]
        B[-1, 0] = 1

        for n in range(Nx - 1):
            A[n, n + 1] = 1
        for n in range(Nx):
            A[-1, n] = -a[Nx - n]
        return cls(A, B, C, D)
        
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

        # Equivalent to _A.sympy.eigenvals()
        
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
    def singular_values(self):
        """List of singular_values."""

        return ExprList(self._A.sympy.singular_values())

    @property    
    def M(self):
        """Modal matrix (eigenvectors of A)."""

        E, L = self._A.diagonalize()
        
        return LaplaceDomainMatrix(E)
    
    @property    
    def Nu(self):
        return self.u.shape[0]    
    
    @property    
    def Nx(self):
        return self.x.shape[0]

    @property    
    def Ny(self):
        return self.y.shape[0]    

    
from .symbols import t, s
from .expr import ExprList
from .sexpr import LaplaceDomainExpression
