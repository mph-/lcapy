"""This module defines the StateSpace class for representing a linear
continuous time-invariant system as a state-space model.

Copyright 2021--2022 Michael Hayes, UCECE

"""

from .sexpr import LaplaceDomainExpression
from .cache import cached_property
from .matrix import Matrix
from .smatrix import LaplaceDomainMatrix
from .tmatrix import TimeDomainMatrix
from .nmatrix import DiscreteTimeDomainMatrix
from .statespacebase import StateSpaceBase
from .texpr import t, texpr
from .expr import expr
from .sym import ssym
import sympy as sym

# Note, need to be careful with simplification since a pole can
# cancel a zero.


class StateSpace(StateSpaceBase):
    """Continuous-time linear time-invariant state space model."""

    @property
    def u(self):
        """Input vector."""
        if self._u is None:
            self._u = TimeDomainMatrix([texpr('u_%d(t)' % i) for i in
                                        range(self.Nu)])
        return self._u

    @property
    def x(self):
        """State variable vector."""
        if self._x is None:
            self._x = TimeDomainMatrix([texpr('x_%d(t)' % i) for i in
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
            self._x0 = TimeDomainMatrix([texpr('x_0_%d(t)' % i) for i in
                                         range(self.Nx)])
        return self._x0

    @property
    def y(self):
        """Output vector."""
        if self._y is None:
            self._y = TimeDomainMatrix([texpr('y_%d(t)' % i) for i in
                                        range(self.Ny)])
        return self._y

    @property
    def U(self):
        """Laplace transform of input vector."""
        return LaplaceDomainMatrix(self.u.LT())

    @property
    def X(self):
        """Laplace transform of state-variable vector."""
        return LaplaceDomainMatrix(self.x.LT())

    @property
    def Y(self):
        """Laplace transform of output vector."""
        return LaplaceDomainMatrix(self.y.LT())

    @cached_property
    def H(self):
        """X(s) / U(s)"""
        return LaplaceDomainMatrix(self.Phi * self._B)

    @property
    def h(self):
        """ILT{X(s) / U(s)}"""
        return TimeDomainMatrix(self.H.ILT(causal=True))

    @cached_property
    def G(self):
        """System transfer functions.
        For a SISO system, use G[0].
        """
        return LaplaceDomainMatrix(self._C * self.H + self._D)

    def state_equations(self):
        """System of first-order differential state equations:

        dotx(t) = A x(t) + B u(t)

        where x is the state vector and u is the input vector.
        """

        return expr(sym.Eq(self.dotx,
                           sym.MatAdd(sym.MatMul(self._A.sympy, self.x.sympy),
                                      sym.MatMul(self._B.sympy, self.u.sympy)),
                           evaluate=False))

    def output_equations(self):
        """System of output equations:

        y(t) = C x(t) + D u(t)

        where y is the output vector, x is the state vector and u is
        the input vector.

        """

        return expr(sym.Eq(self.y,
                           sym.MatAdd(sym.MatMul(self._C.sympy, self.x.sympy),
                                      sym.MatMul(self._D.sympy, self.u.sympy)),
                           evaluate=False))

    @property
    def g(self):
        """System impulse responses."""
        return TimeDomainMatrix(self.G.ILT(causal=True))

    @cached_property
    def Phi(self):
        """s-domain state transition matrix."""

        M = LaplaceDomainMatrix(sym.eye(self.Nx) * ssym - self._A.sympy)
        return LaplaceDomainMatrix(M.inv())

    @cached_property
    def phi(self):
        """State transition matrix."""
        return TimeDomainMatrix(self.Phi.ILT(causal=True))

    def characteristic_polynomial(self):
        """Characteristic polynomial (aka system polynomial).

        `lambda(s) = |s * I - A|`"""

        M = LaplaceDomainMatrix(sym.eye(self.Nx) * ssym - self._A.sympy)
        return LaplaceDomainExpression(M.det()).simplify()

    @cached_property
    def P(self):
        """Characteristic polynomial (aka system polynomial).

        `lambda(s) = |s * I - A|`"""

        return self.characteristic_polynomial()

    @cached_property
    def Lambda(self):
        """Diagonal matrix of eigenvalues."""

        # Perhaps faster to use diagonalize
        # E, L = self.A.diagonalize()
        # return L

        e = self.eigenvalues
        return LaplaceDomainMatrix(sym.diag(*e))

    @cached_property
    def M(self):
        """Modal matrix (eigenvectors of A)."""

        E, L = self._A.sympy.diagonalize()

        return LaplaceDomainMatrix(E)

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

    def generalized_bilinear_transform(self, alpha=0.5):

        from .sym import dt
        from .dtstatespace import DTStateSpace

        if alpha < 0 or alpha > 1:
            raise ValueError("alpha must be between 0 and 1 inclusive")

        I = sym.eye(self.Nx)
        M = I - alpha * dt * self.A
        Minv = M.inv()

        Ad = Minv * (I + (1 - alpha) * dt * self.A)
        Bd = Minv * dt * self.B
        Cd = (Minv.T * self.C.T).T
        Dd = self.D + alpha * self.C * Bd

        # FIXME for u, y, x, x0.
        return DTStateSpace(Ad, Bd, Cd, Dd,
                            DiscreteTimeDomainMatrix(self._u),
                            DiscreteTimeDomainMatrix(self._y),
                            DiscreteTimeDomainMatrix(self._x),
                            DiscreteTimeDomainMatrix(self._x0))

    def discretize(self, method='bilinear', alpha=0.5):
        """Convert to a discrete-time state space approximation.

        The default method is 'bilinear'.  Other methods are
        'forward_euler', 'backward_euler', and 'gbf'.
        The latter has a parameter `alpha`."""

        if method == 'gbf':
            return self.generalized_bilinear_transform(alpha)
        elif method in ('bilinear', 'tustin'):
            return self.generalized_bilinear_transform(0.5)
        elif method in ('euler', 'forward_diff', 'forward_euler'):
            return self.generalized_bilinear_transform(0)
        elif method in ('backward_diff', 'backward_euler'):
            return self.generalized_bilinear_transform(1)
        else:
            raise ValueError('Unsupported method %s' % method)

    @classmethod
    def from_circuit(cls, cct, node_voltages=None, branch_currents=None):

        from .statespacemaker import StateSpaceMaker

        return StateSpaceMaker.from_circuit(cct, node_voltages, branch_currents)
