from .mnacpts import L, C, I, V
from .matrix import Matrix
from .sym import sympify
import sympy as sym

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
        sourcenames = []
        for elt in independent_sources:
            if isinstance(elt, V):
                expr = elt.cpt.voc
            else:
                expr = elt.cpt.isc

            sources.append(expr)
            sourcenames.append(elt.name)

        statesyms = sympify(statenames)
        sourcesyms = sympify(sourcenames)            

        subsdict = {}
        for var, sym1 in zip(statevars, statesyms):
            subsdict[var] = sym1
        for expr, sym1 in zip(sources, sourcesyms):
            subsdict[expr] = sym1            

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
                y.append(sympify('v%s(t)' % node))
                
        Cmat, b = sym.linear_eq_to_matrix(yexprs, *statesyms)
        D, b = sym.linear_eq_to_matrix(yexprs, *sourcesyms)

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

        M = Matrix(sym.eye(len(self.x)) * s - self.A)
        return M.inv()
            
        
from .symbols import t, s
