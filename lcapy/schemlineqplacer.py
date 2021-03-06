from numpy import zeros, hstack, dot, argsort, min, max
from scipy.linalg import lu, inv
from .cnodes import Cnodes
from .schemplacerbase import SchemPlacerBase

class Constraint(object):

    def __init__(self, size, stretch):
        self.size = size
        self.stretch = stretch
        
    def __repr__(self):
        return '%s%s' % (self.size, '*' if self.stretch else '')

        
class SysEq(object):

    def __init__(self, direction, nodes, debug=0):

        self.direction = direction
        self.cnodes = Cnodes(nodes)
        self.debug = debug
        self.constraints = {}
        
    def link(self, n1, n2):
        """Make nodes n1 and n2 share common node"""
        
        self.cnodes.link(n1, n2)        

    def add(self, elt, n1, n2, size, stretch):

        if size == 0:
            return

        if size < 0:
            n1, n2 = n2, n1
            size = -size

        if n1 in self.cnodes:
            n1 = self.cnodes[n1]
        if n2 in self.cnodes:
            n2 = self.cnodes[n2]

        key = n1, n2
        key2 = n2, n1        
        if key not in self.constraints and key2 not in self.constraints:
            self.constraints[key] = Constraint(size, stretch)
            return

        if key2 in self.constraints:
            size = -size
            key = key2
        constraint = Constraint(size, stretch)            
            
        constraint2 = self.constraints[key]
        if not constraint2.stretch:
            if (not constraint.stretch and constraint2.size != constraint.size):
                raise ValueError('Incompatible fixed constraint of size %s and %s' % (constraint.size, constraint2.size))
            self.constraints[key] = constraint
            return

        if constraint.size > constraint2.size:
            self.constraints[key] = constraint            

    def solve(self, stage=None):

        Nstretches = 0
        for key, constraint in self.constraints.items():
            if constraint.stretch:
                Nstretches += 1

        Nnodes = len(self.cnodes)
        Nunknowns = Nnodes - 1 + Nstretches
        Nconstraints = len(self.constraints)

        A = zeros((Nconstraints, Nconstraints + Nnodes - 1))
        b = zeros((Nconstraints, 1))

        cnode_map = {}
        num = 0
        for node, cnode in self.cnodes.items():
            if cnode not in cnode_map:
                cnode_map[cnode] = num
                num += 1

        m = 0
        s = 0

        for key, constraint in self.constraints.items():
            m1 = cnode_map[key[0]]
            m2 = cnode_map[key[1]]
            size = constraint.size

            if m1 != 0:
                A[m, m1 - 1] = -1
            if m2 != 0:        
                A[m, m2 - 1] = 1
            if constraint.stretch:
                A[m, s + Nconstraints - 1] = -1
                s += 1

            b[m] = size
            m += 1

        # Augmented matrix.
        W = hstack((A, b))

        # Extract row-echelon form U.
        PL, U = lu(W, permute_l=True)

        Ur = zeros((U.shape[0], U.shape[0]))

        bound = zeros(A.shape[-1])
        
        c = 0
        for r in range(Ur.shape[0]):
            for c1 in range(c, U.shape[1]):
                if U[r, c1] != 0:
                    Ur[:, r] = U[:, c1]
                    bound[c1] = 1
                    c = c1 + 1
                    break

        br = U[:,-1]
        # A slow one-liner!
        x = dot(inv(Ur), br)[0: Nnodes - 1]

        minx = min(x)
        width = max(x) - minx

        pos = {}
        for node, cnode in self.cnodes.items():
            m = cnode_map[self.cnodes[node]]
            if m == 0:
                pos[node] = 0
            else:
                pos[node] = x[m - 1]                

        self.W = W
        self.b = b

        if ((self.direction == 'horizontal' and not (self.debug & 4))
            or (self.direction == 'vertical' and self.debug & 4)):
        
            if self.debug & 8:
                from .matrix import Matrix
                W = Matrix(W)
                print(W.latex())
                
            if self.debug & 16:
                from .matrix import Matrix
                b = Matrix(b)
                print(b.latex())            

        return pos, width
        

class SchemLineqPlacer(SchemPlacerBase):

    def __init__(self, elements, nodes, debug=0):

        self.elements = elements
        self.nodes = nodes
        self.debug = debug
        
        self.xgraph = SysEq('horizontal', nodes, debug)
        self.ygraph = SysEq('vertical', nodes, debug)

