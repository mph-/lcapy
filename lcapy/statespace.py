"""
This module performs state-space analysis.

Copyright 2019 Michael Hayes, UCECE

"""

from .mnacpts import L, C, I, V
from .matrix import Matrix
from .sym import sympify, ssym
import sympy as sym

# TODO
# 1. Choose better state var names for anonymous C and L. 
# 2. Use a better Matrix class that preserves the class of each
# element, where possible.

# Independent sources
# V1 1 0
# V1 1 0 {10 * u(t)}
# V1 1 0 ac 10
#
# In each case, we use v1(t) as the independent source when determing
# the A, B, C, D matrices.  We can then substitute the known value at
# the end.

class StateSpace(object):
    """This converts a circuit to state-space representation."""

    def __init__(self, cct):

        inductors = []
        capacitors = []
        independent_sources = []

        sscct = cct._new()
        cpt_map = {}

        for key, elt in cct.elements.items():
            ssnet = elt.ss_model()
            sselt = sscct._add(ssnet)
            name = elt.name
            cpt_map[name] = sselt.name
            
            if isinstance(elt, L):
                if sselt.name in cct.elements:
                    raise ValueError('Name conflict %s, either rename the component or iprove the code!' % sselt.name)

                inductors.append(elt)
            elif isinstance(elt, C):
                if sselt.name in cct.elements:
                    raise ValueError('Name conflict %s, either rename the component or iprove the code!' % sselt.name)
                
                capacitors.append(elt)
            elif isinstance(elt, (I, V)):
                independent_sources.append(elt)                
                
        self.cct = cct
        self.sscct = sscct
        
        # Replace inductors with current sources and capacitors with
        # voltage sources.

        dotx_exprs = []
        statevars = []
        statenames = []
        for elt in inductors + capacitors:
            name = cpt_map[elt.name]

            if isinstance(elt, L):
                # Inductors  v = L di/dt  so need v across the L
                expr = sscct[name].v / elt.cpt.L
                var = sscct[name].isc
            else:
                # Capacitors  i = C dv/dt  so need i through the C
                expr = sscct[name].i / elt.cpt.C
                var = sscct[name].voc

            dotx_exprs.append(expr)
            statevars.append(var)
            statenames.append(name)

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

        statesyms = sympify(statenames)
        sourcesyms = sympify(sourcenames)            

        subsdict = {}
        for var, sym1 in zip(statevars, statesyms):
            subsdict[var] = sym1
        for var, sym1 in zip(sourcevars, sourcesyms):
            subsdict[var] = sym1            

        for m, expr in enumerate(dotx_exprs):
            dotx_exprs[m] = expr.subs(subsdict).expr.expand()
                
        A, b = sym.linear_eq_to_matrix(dotx_exprs, *statesyms)
        B, b = sym.linear_eq_to_matrix(dotx_exprs, *sourcesyms)

        # What should be the output vector?  Nodal voltages, branch
        # currents. or both.   Let's start with nodal voltages.

        # Nodal voltages
        ynames = []
        yexprs = []
        y = []

        enodes = list(cct.equipotential_nodes.keys())
        enodes = sorted(enodes)
        for node in enodes:
            if node != '0':
                ynames.append(node)
                yexprs.append(self.sscct.get_vd(node, '0').subs(subsdict).expand())
                y.append(Vt('v%s(t)' % node))

        Cmat, b = sym.linear_eq_to_matrix(yexprs, *statesyms)
        D, b = sym.linear_eq_to_matrix(yexprs, *sourcesyms)

        # Note, Matrix strips the class from each element...
        self.x = Matrix(statevars)

        self.dotx = Matrix([sym.Derivative(x1, t) for x1 in self.x])

        self.u = Matrix(sources)

        self.A = Matrix(A)
        self.B = Matrix(B)        
            
        # Perhaps could use v_R1(t) etc. as the output voltages?
        self.y = Matrix(y)

        self.C = Matrix(Cmat)
        self.D = Matrix(D)

    def state_equations(self):
        """Return system of first-order differential state equations:

        dotx = A x + B u

        where x is the state vector and u is the input vector.
        """
        
        return sym.Eq(self.dotx, sym.MatAdd(sym.MatMul(self.A, self.x), sym.MatMul(self.B, self.u)))

    def output_equations(self):
        """Return system of output equations:

        y = C x + D u

        where y is the output vector, x is the state vector and u is
        the input vector.

        """
        
        return sym.Eq(self.y, sym.MatAdd(sym.MatMul(self.C, self.x), sym.MatMul(self.D, self.u)))

    @property
    def Phi(self):
        """Return s-domain state transition matrix."""

        M = Matrix(sym.eye(len(self.x)) * ssym - self.A)
        return M.inv()

    @property
    def phi(self):
        """Return state transition matrix."""        
        return self.Phi.inverse_laplace(causal=True)
        
    @property
    def U(self):
        """Return Laplace transform of input vector."""
        return self.u.laplace()

    @property
    def X(self):
        """Return Laplace transform of state-variable vector."""        
        return self.x.laplace()

    @property
    def Y(self):
        """Return Laplace transform of output vector."""        
        return self.y.laplace()    

    @property
    def H(self):
        """X(s) / U(s)"""

        return self.Phi * self.B

    @property
    def G(self):
        """Return system transfer functions: Y(s) / U(s)"""

        return self.C * self.H + self.D

    def characteristic_polynomial(self):
        """Return characteristic polynomial (aka system polynomial).

        lambda(s) = |s * I - A|
        """

        M = Matrix(sym.eye(len(self.x)) * ssym - self.A)        
        return sExpr(M.det()).simplify()

    @property
    def P(self):
        """Return characteristic polynomial (aka system polynomial).

        lambda(s) = |s * I - A|
        """        
        return self.characteristic_polynomial()

    @property        
    def eigenvalues_dict(self):
        """Return dictionary of eigenvalues, the roots of the characteristic
        polynomial (equivalent to the poles of Phi(s)).  The
        dictionary values are the multiplicity of the eigenvalues.

        For a list of eigenvalues use eigenvalues.

        """        

        return self.characteristic_polynomial().roots()
        
    @property        
    def eigenvalues(self):
        """Return list of eigenvalues, the roots of the characteristic polynomial
        (equivalent to the poles of Phi(s))."""
        
        roots = self.eigenvalues_dict
        e = []

        # Replicate duplicated eigenvalues and return as a list.
        for v, n in roots.items():
            for m in range(n):
                e.append(v)
        return e

    @property    
    def Lambda(self):
        """Diagonal matrix of eigenvalues."""

        # Perhaps faster to use diagonalize
        # E, L = self.A.diagonalize()
        # return L
        
        e = self.eigenvalues
        
        M = Matrix(sym.diag(*e))
        return  M

    @property        
    def eigenvectors(self):
        """Return list of tuples (eigenvalue, multiplicity of igenvalue,
        basis of the eigenspace) of A.

        """
        return self.A.eigenvects()
    
    @property    
    def M(self):
        """Return modal matrix (eigenvectors of A)."""

        E, L = self.A.diagonalize()
        
        return Matrix(E)
    
    
from .symbols import t, s
from .texpr import Vt, tExpr
from .sexpr import sExpr
