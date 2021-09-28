"""This module defines the StateSpace class for representing a linear
time-invariant system as a state-space model.

Copyright 2019-2021 Michael Hayes, UCECE

"""

from .matrix import Matrix
from .vector import Vector
from .smatrix import LaplaceDomainMatrix
from .tmatrix import TimeDomainMatrix
from .sym import ssym
from .texpr import t, texpr
from .expr import expr
from .cache import cached_property
import sympy as sym


class StateSpace(object):

    def __init__(self, A, B, C, D, u=None, y=None, x=None, x0=None):
        """Create continuous state-space object where:

        A is Nx x Nx state matrix
        B is Nx x Nu input matrix
        C is Ny x Nx output matrix
        D is Ny x Nu feedthrough matrix

        u is Nu x 1 input vector
        x is Nx x 1 state vector
        x0 is Nx x 1 state initial vector
        y is Ny x 1 state output vector
        """

        A = Matrix(A)
        B = Matrix(B)
        C = Matrix(C)
        D = Matrix(D)

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

        if x is not None:
            if x.shape[0] != Nx:
                raise ValueError('x vector has wrong dimension')
            x = Vector(x)            
        if x0 is not None:
            if x0.shape[0] != Nx:
                raise ValueError('x0 vector has wrong dimension')
            x0 = Vector(x0)                        
        if u is not None:
            if u.shape[0] != Nu:
                raise ValueError('u vector has wrong dimension')
            u = Vector(u)
        if y is not None:
            if y.shape[0] != Ny:
                raise ValueError('y vector has wrong dimension')
            y = Vector(y)
            
        self._A = A
        self._B = B
        self._C = C
        self._D = D

        self._u = u
        self._x = x
        self._x0 = x0        
        self._y = y

    @classmethod
    def from_ba(cls, b, a, form='CCF'):
        """Create state-space representation from transfer function
        specified with numerator and denominator coefficients.

         Note, state-space representations are not unique and are
        determined by the `form` argument.  Currently this can be
        'CCF' for the controllable canonical form.

        """
        return cls.from_transfer_function_coeffs(cls, b, a, form)
        
    @classmethod
    def from_transfer_function_coeffs(cls, b, a, form='CCF'):    
        """Create state-space representation from transfer function
        specified with numerator and denominator coefficients.

         Note, state-space representations are not unique and are
        determined by the `form` argument.  Currently this can be
        'CCF' for the controllable canonical form.

        """        

        if form == 'CCF':
            return cls.from_ba_CCF(b, a)
        elif form == 'OCF':
            return cls.from_ba_OCF(b, a)
        elif form == 'DCF':
            return cls.from_ba_DCF(b, a)        
        else:
            raise ValueError('Only CCF, DCF, and OCF forms are currently supported')

    @classmethod
    def from_ba_CCF(cls, b, a):
        
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
            # Need extended state-space representation...
            raise ValueError('Improper transfer function; require derivatives of input')

        Nx = len(a) - 1
        Nu = 1
        Ny = 1

        A = Matrix.zeros(Nx, Nx)
        B = Matrix.zeros(Nx, 1)
        C = Matrix.zeros(1, Nx)
        D = Matrix.zeros(1, 1)

        D[0, 0] = b[0]
        for n in range(Nx):
            C[0, n] = b[Nx - n] - a[Nx - n] * b[0]
        B[-1, 0] = sym.S.One

        for n in range(Nx - 1):
            A[n, n + 1] = sym.S.One
        for n in range(Nx):
            A[-1, n] = -a[Nx - n]
        return cls(A, B, C, D)

    @classmethod
    def from_ba_OCF(cls, b, a):

        # Aobs = Acon.T
        # Bobs = Ccon.T
        # Cobs = Bcon.T
        # Dobs = Dcon
        
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
            # Need extended state-space representation...
            raise ValueError('Improper transfer function; require derivatives of input')

        Nx = len(a) - 1
        Nu = 1
        Ny = 1

        A = Matrix.zeros(Nx, Nx)
        B = Matrix.zeros(Nx, 1)
        C = Matrix.zeros(1, Nx)
        D = Matrix.zeros(1, 1)

        D[0, 0] = b[0]
        for n in range(Nx):
            B[n, 0] = b[n + 1] - a[n + 1] * b[0]
        C[0, 0] = sym.S.One

        for n in range(Nx - 1):
            A[n, n + 1] = sym.S.One
        for n in range(Nx):
            A[n, 0] = -a[n + 1]
        return cls(A, B, C, D)

    @classmethod
    def from_ba_DCF(cls, b, a):

        from .sexpr import tf
        
        b = list(b)
        a = list(a)
        Nb = len(b)
        Na = len(a)

        if Nb > Na:
            # Need extended state-space representation...
            raise ValueError('Improper transfer function; require derivatives of input')
        
        H = tf(b, a)

        poles = H._ratfun.poles()
        for p in poles:
            if p.n != 1:
                raise ValueError('Require unique poles')
        
        Nx = len(a) - 1
        Nu = 1
        Ny = 1

        A = Matrix.zeros(Nx, Nx)
        B = Matrix.ones(Nx, 1)
        C = Matrix.zeros(1, Nx)
        D = Matrix.zeros(1, 1)

        # FIXME
        if Na == Nb:
            D[0, 0] = b[0]
        else:
            D[0, 0] = 0            

        for n, p in enumerate(poles):
            A[n, n] = p.expr
            C[n] = H._ratfun.residue(p.expr, poles)
        
        return cls(A, B, C, D)        

    def state_equations(self):
        """System of first-order differential state equations:

        dotx = A x + B u

        where x is the state vector and u is the input vector.
        """
        
        return expr(sym.Eq(self.dotx, sym.MatAdd(sym.MatMul(self._A, self.x),
                                                 sym.MatMul(self._B, self.u)),
                           evaluate=False))

    def output_equations(self):
        """System of output equations:

        y = C x + D u

        where y is the output vector, x is the state vector and u is
        the input vector.

        """
        
        return expr(sym.Eq(self.y, sym.MatAdd(sym.MatMul(self._C, self.x),
                                              sym.MatMul(self._D, self.u)),
                           evaluate=False))

    @property
    def u(self):
        """Input vector."""
        if self._u is None:
            self._u = TimeDomainMatrix([texpr('u_%d(t)' % n) for n in
                                        range(self.Nu)])
        return self._u
    
    @property
    def x(self):
        """State variable vector."""
        if self._x is None:
            self._x = TimeDomainMatrix([texpr('x_%d(t)' % n) for n in
                                        range(self.Nx)])
        return self._x

    @cached_property
    def dotx(self):
        """Time derivative of state variable vector."""
        return TimeDomainMatrix([sym.Derivative(x1, t) for x1 in self.x])
    
    @property
    def x0(self):
        """State variable initial value vector."""
        if self._x0 is None:
            self._x0 = TimeDomainMatrix([texpr('x_0_%d(t)' % n) for n in
                                         range(self.Nx)])
        return self._x0        

    @property
    def y(self):
        """Output vector."""
        if self._y is None:
            self._y = TimeDomainMatrix([texpr('y_%d(t)' % n) for n in
                                        range(self.Ny)])
        return self._y    

    @property
    def A(self):
        """State matrix."""
        return self._A

    @property
    def state_matrix(self):
        """State matrix."""
        return self._A    

    @property
    def B(self):
        """Input matrix."""
        return self._B

    @property
    def input_matrix(self):
        """Input matrix."""
        return self._B
    
    @property
    def C(self):
        """Output matrix."""
        return self._C

    @property
    def output_matrix(self):
        """Output matrix."""
        return self._C    

    @property
    def D(self):
        """Feed-through matrix."""
        return self._D

    @property
    def feedthrough_matrix(self):
        """Feed-through matrix."""
        return self._D            

    @cached_property
    def Phi(self):
        """s-domain state transition matrix."""

        M = LaplaceDomainMatrix(sym.eye(self.Nx) * ssym - self._A)
        return LaplaceDomainMatrix(M.inv().canonical())

    @cached_property
    def phi(self):
        """State transition matrix."""        
        return TimeDomainMatrix(self.Phi.inverse_laplace(causal=True))

    @property
    def state_transition_matrix(self):
        """State transition matrix."""        
        return self.phi
        
    @property
    def U(self):
        """Laplace transform of input vector."""
        return LaplaceDomainMatrix(self.u.laplace())

    @property
    def X(self):
        """Laplace transform of state-variable vector."""        
        return LaplaceDomainMatrix(self.x.laplace())

    @property
    def Y(self):
        """Laplace transform of output vector."""        
        return LaplaceDomainMatrix(self.y.laplace())

    @cached_property
    def H(self):
        """X(s) / U(s)"""
        return LaplaceDomainMatrix(self.Phi * self._B).canonical()

    @property
    def h(self):
        """ILT{X(s) / U(s)}"""        
        return TimeDomainMatrix(self.H.inverse_laplace(causal=True))

    @cached_property
    def G(self):
        """System transfer functions.
        For a SISO system, use G[0].
        """
        return LaplaceDomainMatrix(self._C * self.H + self._D).canonical()

    @property
    def transfer_functions(self):
        """System transfer functions.  See also `transfer_function`"""
        return self.G
        
    @property
    def transfer_function(self):
        """System transfer function for a SISO system.  See also
        `transfer_functions`"""
        return self.G[0]
    
    @property
    def g(self):
        """System impulse responses."""        
        return TimeDomainMatrix(self.G.inverse_laplace(causal=True))
    
    def characteristic_polynomial(self):
        """Characteristic polynomial (aka system polynomial).

        `lambda(s) = |s * I - A|`"""

        M = Matrix(sym.eye(self.Nx) * ssym - self._A)        
        return LaplaceDomainExpression(M.det()).simplify()

    @cached_property
    def P(self):
        """Characteristic polynomial (aka system polynomial).

        `lambda(s) = |s * I - A|`"""        

        return self.characteristic_polynomial().canonical()

    @cached_property        
    def eigenvalues_dict(self):
        """Dictionary of eigenvalues, the roots of the characteristic
        polynomial (equivalent to the poles of Phi(s)).  The
        dictionary values are the multiplicity of the eigenvalues.

        For a list of eigenvalues use eigenvalues."""        

        # Equivalent to _A.sympy.eigenvals()
        
        return self.characteristic_polynomial().roots()
        
    @cached_property        
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

    @cached_property    
    def Lambda(self):
        """Diagonal matrix of eigenvalues."""

        # Perhaps faster to use diagonalize
        # E, L = self.A.diagonalize()
        # return L
        
        e = self.eigenvalues
        return LaplaceDomainMatrix(sym.diag(*e))

    @cached_property        
    def eigenvectors(self):
        """List of tuples (eigenvalue, multiplicity of eigenvalue,
        basis of the eigenspace) of A."""
        
        return self._A.eigenvects()
    
    @cached_property        
    def singular_values(self):
        """List of singular_values."""

        return ExprList(self._A.sympy.singular_values())

    @cached_property    
    def M(self):
        """Modal matrix (eigenvectors of A)."""

        E, L = self._A.diagonalize()
        
        return LaplaceDomainMatrix(E)
    
    @property
    def Nx(self):
        """Number of state variables (the system order)."""
        return self._A.shape[0]

    @property
    def Nu(self):
        """Number of inputs."""    
        return self._B.shape[1]

    @property
    def Ny(self):
        """Number of outputs."""        
        return self._C.shape[0]
    
    @cached_property    
    def is_symbolic(self):

        return ((self.A.symbols != {}) or (self.B.symbols != {}) or
                (self.C.symbols != {}) or (self.D.symbols != {}))
    
    @cached_property    
    def is_stable(self):
        """True if system is stable."""

        if self.is_symbolic:
            return None
        
        for e in self.eigenvalues:
            if e.real > 0:
                return False
        return True

    @cached_property        
    def controllability_matrix(self):
        """Return controllability matrix."""
        
        B = self.B
        A = self.A

        R = B
        Q = B
        for n in range(self.Nx - 1):
            Q = A * Q
            R = R.hstack(R, Q)
        return R

    @property        
    def R(self):
        """Return controllability matrix."""

        return self.controllability_matrix    

    @cached_property        
    def is_controllable(self):
        """True if system is controllable."""        

        R = self.controllability_matrix
        return R.rank() == R.shape[0]

    @cached_property        
    def observability_matrix(self):
        """Return observability matrix."""        

        C = self.C
        A = self.A

        O = C
        Q = C
        for n in range(self.Nx - 1):
            Q *= A
            O = O.vstack(O, Q)
        return O    

    @property        
    def O(self):
        """Return observability matrix."""

        return self.observability_matrix
    
    @cached_property        
    def is_observable(self):
        """True if system is observable."""                

        O = self.observability_matrix
        return O.rank() == O.shape[1]

    @cached_property            
    def controllability_gramian(self):
        """Controllability gramian matrix."""

        from scipy import linalg

        B = self.B.evaluate()
        Q = -B @ B.T

        # Find Wc given A @ Wc + Wc @ A.T = Q        
        # Wc > o if (A, B) controllable
        Wc = linalg.solve_continuous_lyapunov(self.A.evaluate(), Q)

        # Wc should be symmetric positive semi-definite
        Wc = (Wc + Wc.T) / 2        
        
        return Matrix(Wc)

    @property            
    def Wc(self):
        """Controllability gramian matrix."""
        return self.controllability_gramian
    
    @property            
    def reachability_gramian(self):
        """Reachability gramian matrix.  This is equivalent to the
        controllability gramian matrix for a linear time independent
        system."""

        return self.controllability_gramian

    @property            
    def Wr(self):
        """Reachability gramian matrix."""
        return self.reachability_gramian
    
    @cached_property            
    def observability_gramian(self):
        """Observability gramian matrix."""

        from scipy import linalg

        C = self.C.evaluate()
        Q = -C.T @ C

        # Find Wo given A.T @ Wo + Wo @ A = Q
        # Wo > o if (C, A) observable
        Wo = linalg.solve_continuous_lyapunov(self.A.evaluate().T, Q)

        # Wo should be symmetric positive semi-definite
        Wo = (Wo + Wo.T) / 2
        
        return Matrix(Wo)

    @property            
    def Wo(self):
        """Observability gramian matrix."""
        return self.observability_gramian
    
    @cached_property            
    def hankel_singular_values(self):

        from numpy import sqrt
        from numpy.linalg import eig
        
        Wc = self.controllability_gramian.evaluate()
        Wo = self.observability_gramian.evaluate()

        e, UT = eig(Wc @ Wo)

        return expr(sqrt(e), rational=False)

    @cached_property                
    def balanced_transformation(self):
        """Return the transformation matrix `T` required to balance the
        controllability and observability gramians.
        
        `Wob = Tinv.T * Wo * Tinv
        Wcb = T * Wc * T.T`

        where `Tinv = T.inv()`

        """

        from scipy import linalg
        from numpy import sqrt, diag

        if not self.is_stable:
            raise ValueError('System not stable')
        
        Wc = self.controllability_gramian.evaluate()
        Wo = self.observability_gramian.evaluate()

        # Wc = R.T @ R
        R = linalg.cholesky(Wc)

        Y = R @ Wo @ R.T
        # Y = U @ diag(e) @ U.T
        e, U = linalg.eig(Y, left=True, right=False)

        # e is a vector of squared Hankel singular values
        Einv = diag(e ** -0.25)
        
        Tinv = R.T @ U @ Einv
        return Matrix(Tinv).inv()

    def balance(self):
        """Return new StateSpace object that has the controllability and
        observability gramians equal to the diagonal matrix with the
        Hankel singular values on the diagonal."""

        T = self.balanced_transformation
        return self.transform(T)
    
    def transform(self, T):

        Tinv = T.inv()
        Ap = T * self.A * Tinv
        Bp = T * self.B
        Cp = self.C * Tinv
        
        return self.__class__(Ap, Bp, Cp, self.D,
                              self._u, self._y, self._x, self._x0)

    def reduce(self, elim_states=None, method='truncate'):
        """Perform model reduction given array `elim_states` of
        states to remove."""

        from numpy import array, sort, hstack, linalg

        # Perhaps also allow a Boolean array to select states.
        
        melim = sort(elim_states)
        mkeep = []

        for i in range(0, self.Nx):
            if i not in melim:
                mkeep.append(i)

        # A1 is a matrix of all columns of A to keep
        A1 = self.A[:, mkeep[0]]
        for i in mkeep[1:]:
            A1 = hstack((A1, self.A[:, i]))
        A11 = A1[mkeep, :]
        A21 = A1[melim, :]
        
        # A2 is a matrix of all columns of A to eliminate
        A2 = self.A[:, melim[0]]
        for i in melim[1:]:
            A2 = hstack((A2, self.A[:, i]))
        A12 = A2[mkeep, :]
        A22 = A2[melim, :]
        
        C1 = self.C[:, mkeep]
        C2 = self.C[:, melim]
        B1 = self.B[mkeep, :]
        B2 = self.B[melim, :]
        D = self.D
        
        if method == 'truncate':
            Ar = A11
            Br = B1
            Cr = C1
            Dr = D
        elif method == 'matchdc':
            A22I = linalg.inv(A22)
            
            Ar = A11 - A12 * A22I * A21
            Br = B1 - A12 * A22I * B2
            Cr = C1 - C2 * A22I * A21
            Dr = D - C2 * A22I * B2
        else:
            raise ValueError("Reduction method %s is not supported.  Try 'matchdc' or 'truncate'" % method)

        x = self._x
        x0 = self._x0
        u = self._u
        y = self._y
        if x is not None:
            x = array(x)[mkeep]
        if x0 is not None:
            x0 = array(x0)[mkeep]
        if u is not None:
            u = array(u)[mkeep]
        if y is not None:
            y = array(y)[mkeep]            
        
        return self.__class__(Ar, Br, Cr, Dr, u, y, x, x0)        

    def balance_reduce(self, threshold, method='truncate'):
        """Perform balanced model reduction where the states with hankel
        singular values smaller than `threshold` are removed."""

        from numpy import arange
        
        hsv = self.hankel_singular_values

        elim_states = arange(self.Nx)[hsv.numpy.squeeze() < threshold]
        return self.reduce(elim_states, method)
    
from .symbols import t, s
from .expr import ExprList
from .sexpr import LaplaceDomainExpression
