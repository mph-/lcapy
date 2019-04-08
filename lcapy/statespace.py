from .mnacpts import L, C, I, V
from .matrix import Matrix
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

        statesyms = sym.sympify(statenames)
        sourcesyms = sym.sympify(sourcenames)            

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
        for node in cct.equipotential_nodes.keys():
            if node != '0':
                ynames.append(node)
                yexprs.append(self.sscct.get_vd(node, '0').subs(subsdict).expand())

        Cmat, b = sym.linear_eq_to_matrix(yexprs, *statesyms)
        D, b = sym.linear_eq_to_matrix(yexprs, *sourcesyms)                

        self.C = Matrix(Cmat)
        self.D = Matrix(D)

        self.x = Matrix(statevars)        
        self.u = Matrix(sources)

        self.dotx = Matrix([sym.Derivative(sym1, t) for sym1 in statevars])
    
from .symbols import t
