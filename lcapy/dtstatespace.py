"""This module defines the StateSpace class for representing a linear
discrete time-invariant system as a state-space model.

Copyright 2021 Michael Hayes, UCECE

"""

from .cache import cached_property
from .matrix import Matrix
from .zmatrix import ZDomainMatrix
from .nmatrix import DiscreteTimeDomainMatrix
from .statespacebase import StateSpaceBase
from .nexpr import n, nexpr
from .expr import expr
from .sym import zsym
import sympy as sym

# TODO: reachability implies controllability but controllability only
# implies reachability if A matrix is full rank.


class DTStateSpace(StateSpaceBase):
    """Discrete-time linear time-invariant state space model."""

    @property
    def u(self):
        """Input vector."""
        if self._u is None:
            self._u = DiscreteTimeDomainMatrix([nexpr('u_%d(n)' % i) for i in
                                                range(self.Nu)])
        return self._u

    @property
    def x(self):
        """State variable vector."""
        if self._x is None:
            self._x = DiscreteTimeDomainMatrix([nexpr('x_%d(n)' % i) for i in
                                                range(self.Nx)])
        return self._x

    @cached_property
    def xnext(self):
        """Time derivative of state variable vector."""
        return DiscreteTimeDomainMatrix([nexpr('x_%d(n + 1)' % i) for i in
                                         range(self.Nx)])

    @property
    def x0(self):
        """State variable initial value vector."""
        if self._x0 is None:
            self._x0 = DiscreteTimeDomainMatrix([nexpr('x_0_%d(n)' % i) for i in
                                                 range(self.Nx)])
        return self._x0

    @property
    def y(self):
        """Output vector."""
        if self._y is None:
            self._y = DiscreteTimeDomainMatrix([nexpr('y_%d(n)' % i) for i in
                                                range(self.Ny)])
        return self._y

    @property
    def U(self):
        """Z transform of input vector."""
        return ZDomainMatrix(self.u.ZT())

    @property
    def X(self):
        """Z transform of state-variable vector."""
        return ZDomainMatrix(self.x.ZT())

    @property
    def Y(self):
        """Z transform of output vector."""
        return ZDomainMatrix(self.y.ZT())

    @cached_property
    def H(self):
        """X(z) / U(z)"""
        return ZDomainMatrix(self.Phi * self._B)

    @property
    def h(self):
        """ILT{X(z) / U(z)}"""
        return DiscreteTimeDomainMatrix(self.H.ILT(causal=True))

    @cached_property
    def G(self):
        """System transfer functions.
        For a SISO system, use G[0].
        """
        return ZDomainMatrix(self._C * self.H + self._D)

    def state_equations(self):
        """System of first-order differential state equations:

        x[n + 1] = A x[n] + B u[n]

        where x is the state vector and u is the input vector.
        """

        return expr(sym.Eq(self.xnext, sym.MatAdd(sym.MatMul(self._A, self.x),
                                                  sym.MatMul(self._B, self.u)),
                           evaluate=False))

    def output_equations(self):
        """System of output equations:

        y[n] = C x[n] + D u[n]

        where y is the output vector, x is the state vector and u is
        the input vector.

        """

        return expr(sym.Eq(self.y, sym.MatAdd(sym.MatMul(self._C, self.x),
                                              sym.MatMul(self._D, self.u)),
                           evaluate=False))

    @property
    def g(self):
        """System impulse responses."""
        return DiscreteTimeDomainMatrix(self.G.IZT(causal=True))

    @cached_property
    def Phi(self):
        """z-domain state transition matrix."""

        M = ZDomainMatrix(sym.eye(self.Nx) * zsym - self._A.sympy)
        return ZDomainMatrix(M.inv())

    @cached_property
    def phi(self):
        """State transition matrix."""
        return DiscreteTimeDomainMatrix(self.Phi.ILT(causal=True))

    def characteristic_polynomial(self):
        """Characteristic polynomial (aka system polynomial).

        `lambda(z) = |z * I - A|`"""

        M = ZDomainMatrix(sym.eye(self.Nx) * zsym - self._A.sympy)
        return ZDomainExpression(M.det()).simplify()

    @cached_property
    def P(self):
        """Characteristic polynomial (aka system polynomial).

        `lambda(z) = |z * I - A|`"""

        return self.characteristic_polynomial()

    @cached_property
    def Lambda(self):
        """Diagonal matrix of eigenvalues."""

        # Perhaps faster to use diagonalize
        # E, L = self.A.diagonalize()
        # return L

        e = self.eigenvalues
        return ZDomainMatrix(sym.diag(*e))

    @cached_property
    def M(self):
        """Modal matrix (eigenvectors of A)."""

        E, L = self._A.diagonalize()

        return ZDomainMatrix(E)

    @cached_property
    def controllability_gramian(self):
        """Controllability gramian matrix."""

        from scipy import linalg

        B = self.B.evaluate()
        Q = -B @ B.T

        # Find Wc given A @ Wc + Wc @ A.T = Q
        # Wc > o if (A, B) controllable
        Wc = linalg.solve_discrete_lyapunov(self.A.evaluate(), Q)

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
        system provided A is not singular."""

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
        Wo = linalg.solve_discrete_lyapunov(self.A.evaluate().T, Q)

        # Wo should be symmetric positive semi-definite
        Wo = (Wo + Wo.T) / 2

        return Matrix(Wo)

    @property
    def Wo(self):
        """Observability gramian matrix."""
        return self.observability_gramian

    def minimum_energy_input(self, steps, xfinal, xinitial=None):
        """Determine minimum energy (least norm) input vector that results in
        state `xfinal` in `steps` time steps from initial state
        `xinitial`."""

        if xinitial is None:
            xinitial = self.x0
        else:
            xinitial = Matrix(xinitial)

        xfinal = Matrix(xfinal)

        if steps < 0:
            raise ValueError('steps must be positive')

        C = self.controllability_matrix_steps(steps)

        xdiff = xfinal - xinitial

        uflip = C.T * (C * C.T).inv() * xdiff
        return uflip[::-1, :]

    def minimum_energy(self, steps, xfinal, xinitial=None):
        """Determine minimum energy that results in state `xfinal` in `steps`
        time steps from initial state `xinitial`."""

        if xinitial is None:
            xinitial = self.x0
        else:
            xinitial = Matrix(xinitial)

        xfinal = Matrix(xfinal)

        if steps < 0:
            raise ValueError('steps must be positive')

        C = self.controllability_matrix_steps(steps)

        xdiff = xfinal - xinitial

        return (xdiff.T * (C * C.T).inv() * xdiff)[0]

    def state_transfer(self, u, xinitial=None):
        """Return transitioned state given inputs `u` specified as a list of
        column vectors; one for each time step."""

        if xinitial is None:
            xinitial = self.x0
        else:
            xinitial = Matrix(xinitial)

        u = Matrix(u)
        steps = u.shape[0]

        C = self.controllability_matrix_steps(steps)

        xfinal = self.A**steps * xinitial + C * u[::-1, :]
        return xfinal


from .zexpr import ZDomainExpression  # nopep8
