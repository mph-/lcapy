from .mnacpts import L, C, I, V
from .matrix import Matrix
from .sym import sympify
import sympy as sym

class StateSpace(object):

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
                inductors.append(elt)
            elif isinstance(elt, C):
                capacitors.append(elt)
            elif isinstance(elt, (I, V)):
                independent_sources.append(elt)                
                
        self.inductors = inductors
        self.capacitors = capacitors
        self.independent_sources = independent_sources

        # Replace inductors with current sources and capacitors with
        # voltage sources.
        self.sscct = sscct

        # Capacitors  i = C dv/dt  so need i through the C
        # Inductors  v = L di/dt  so need v across the L
        dotx_exprs = []
        statevars = []
        statenames = []
        for elt in inductors + capacitors:
            name = cpt_map[elt.name]

            if isinstance(elt, L):
                expr = sscct[name].v / elt.cpt.L
                var = sscct[name].isc
            else:
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

        self.A = Matrix(A)
        self.B = Matrix(B)        
            
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

        self.y = Matrix(y)

        self.C = Matrix(Cmat)
        self.D = Matrix(D)

        self.x = Matrix(statevars)        
        self.u = Matrix(sources)

        self.dotx = Matrix([sym.Derivative(sym1, t) for sym1 in statevars])

    @property
    def state_equations(self):
        """Return system of first-order differential state equations.

        dotx = A x + B u
        """
        
        return sym.Eq(self.dotx, self.A * self.x + self.B * self.u)

    @property
    def output_equations(self):
        """Return system of output equations.

        y = C x + Du
        """
        
        return sym.Eq(self.y, self.C * self.x + self.D * self.u)    

    @property
    def Phi(self):
        """Return s-domain state transition matrix."""

        M = Matrix(sym.eye(len(self.x)) * s - self.A)
        return M.inv()
            
        
from .symbols import t, s
