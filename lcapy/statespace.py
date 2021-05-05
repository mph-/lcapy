"""
This module performs state-space analysis.

Copyright 2019-2020 Michael Hayes, UCECE

"""

from .mnacpts import I, V, C, L
from .expr import expr
from .matrix import Matrix
from .smatrix import LaplaceDomainMatrix
from .tmatrix import TimeDomainMatrix
from .voltage import voltage
from .current import current
from .sym import sympify, ssym
import sympy as sym

__all__ = ('StateSpace', )

# TODO
# 1. Use a better Matrix class that preserves the class of each
# element, where possible.
# 2. Handle namespaces

# Independent sources
# V1 1 0
# V1 1 0 {10 * u(t)}
# V1 1 0 ac 10
#
# In each case, we use v1(t) as the independent source when determining
# the A, B, C, D matrices.  We can then substitute the known value at
# the end.

# There is a generalised (but less common) state-space representation:
#
# dx/dt = A x + B u + C du/dt
# y = D x + E u + F du/dt
#
# This is required when have a current source in series with an inductor;
# here the inductor current is not a state variable since it is equivalent
# to an input.  However, the output equation for the inductor voltage
# requires the derivative of the input current.


def _hack_vars(exprs):
    """Substitute i_Canon1(t) with i_C(t) etc. provided
    there is no i_Canon2(t)."""

    for m, expr1 in enumerate(exprs):
        for c in ('i_V', 'i_C', 'i_L', 'v_C'):
            sym1 = sympify(c + 'anon1(t)')
            sym2 = sympify(c + 'anon2(t)')            
            if expr1.has(sym1) and not expr1.has(sym2):
                expr1 = expr1.subs(sym1, sympify(c + '(t)'))
                exprs[m] = expr1

class StateSpace(object):
    """This converts a circuit to state-space representation.

    It can take a long time for a symbolic circuit with many reactive
    components.

    The currents through inductors and the voltage across capacitors
    are chosen as the state variables.

    This does not (yet) look for degenerate circuits.  These are
    circuits with a loop consisting only of voltage sources and/or
    capacitors, or a cut set consisting only of current sources and/or
    inductors.  One hack is to call simplify() first to remove series
    inductors and parallel capacitors.

    """

    def __init__(self, cct, node_voltages=True, branch_currents=False):

        if not node_voltages and not branch_currents:
            raise ValueError('No outputs')
        
        inductors = []
        capacitors = []
        independent_sources = []

        # Determine state variables (current through inductors and
        # voltage across acapacitors) and replace inductors with
        # current sources and capacitors with voltage sources.
        sscct = cct._new()
        cpt_map = {}

        for key, elt in cct.elements.items():
            ssnet = elt._ss_model()
            sselt = sscct._add(ssnet)
            name = elt.name
            cpt_map[name] = sselt.name
            
            if elt.is_inductor:
                if sselt.name in cct.elements:
                    raise ValueError('Name conflict %s, either rename the component or improve the code!' % sselt.name)

                inductors.append(elt)
            elif elt.is_capacitor:
                if sselt.name in cct.elements:
                    raise ValueError('Name conflict %s, either rename the component or improve the code!' % sselt.name)
                
                capacitors.append(elt)
            elif isinstance(elt, (I, V)):
                independent_sources.append(elt)                
                
        self.cct = cct
        self.sscct = sscct
        # sscct can be analysed in the time domain since it has no
        # reactive components.  However, for large circuits
        # this can take a long time due to inversion of the MNA matrix.
        
        dotx_exprs = []
        statevars = []
        statenames = []
        initialvalues = []
        for elt in inductors + capacitors:
            name = cpt_map[elt.name]

            if isinstance(elt, L):
                # Inductors  v = L di/dt  so need v across the L
                expr = sscct[name].v / elt.cpt.L
                var = -sscct[name].isc
                x0 = elt.cpt.i0
            else:
                # Capacitors  i = C dv/dt  so need i through the C
                # The current is negated since it is from a source V_Cx
                expr = -sscct[name].i / elt.cpt.C
                var = sscct[name].voc
                x0 = elt.cpt.v0

            dotx_exprs.append(expr)
            statevars.append(var)
            statenames.append(name)
            initialvalues.append(x0)

        statesyms = sympify(statenames)

        # Determine independent sources.
        sources = []
        sourcevars = []
        sourcenames = []
        for elt in independent_sources:
            name = cpt_map[elt.name]
            
            if isinstance(elt, V):
                expr = elt.cpt.voc
                var = sscct[name].voc                
            else:
                expr = elt.cpt.isc
                var = sscct[name].isc

            sources.append(expr)
            sourcevars.append(var)
            sourcenames.append(name)

        sourcesyms = sympify(sourcenames)            

        # linear_eq_to_matrix expects only Symbols and not AppliedUndefs,
        # so substitute.
        subsdict = {}
        for var, sym1 in zip(statevars, statesyms):
            subsdict[var.expr] = sym1
        for var, sym1 in zip(sourcevars, sourcesyms):
            subsdict[var.expr] = sym1       

        for m, expr in enumerate(dotx_exprs):
            dotx_exprs[m] = expr.subs(subsdict).expr.expand()

        A, b = sym.linear_eq_to_matrix(dotx_exprs, *statesyms)
        if sourcesyms != []:
            B, b = sym.linear_eq_to_matrix(dotx_exprs, *sourcesyms)
        else:
            B = []

        # Determine output variables.
        yexprs = []
        y = []

        if node_voltages:
            for node in cct.node_list:
                if node == '0':
                    continue
                yexprs.append(self.sscct[node].v.subs(subsdict).expand().expr)
                # Note, this can introduce a name conflict
                y.append(voltage('v_%s(t)' % node))

        if branch_currents:
            for name in cct.branch_list:
                # Perhaps ignore L since the current through it is a
                # state variable?
                name2 = cpt_map[name]                    
                yexprs.append(self.sscct[name2].i.subs(subsdict).expand().expr)
                y.append(current('i_%s(t)' % name))                    

        Cmat, b = sym.linear_eq_to_matrix(yexprs, *statesyms)
        if sourcesyms != []:        
            D, b = sym.linear_eq_to_matrix(yexprs, *sourcesyms)
        else:
            D = []

        # Rewrite vCanon1(t) as vC(t) etc if appropriate.
        _hack_vars(statevars)
        _hack_vars(sources)
        
        # Note, Matrix strips the class from each element...
        self._x = TimeDomainMatrix(statevars)

        self._x0 = Matrix(initialvalues)
        
        self.dotx = TimeDomainMatrix([sym.Derivative(x1, t) for x1 in self._x])

        self._u = TimeDomainMatrix(sources)

        self._A = Matrix(A)
        self._B = Matrix(B)        
            
        # Perhaps could use v_R1(t) etc. as the output voltages?
        self._y = TimeDomainMatrix(y)

        self._C = Matrix(Cmat)
        self._D = Matrix(D)

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
