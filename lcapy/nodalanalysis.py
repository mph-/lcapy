"""This module performs nodal analysis.  It is primarily for showing
the equations rather than for evaluating them.

Copyright 2019-202 Michael Hayes, UCECE

"""

from .circuitgraph import CircuitGraph
from .tmatrix import tMatrix
from .utils import equation
import sympy as sym

__all__ = ('NodalAnalysis', )


class NodalAnalysis(object):
    """
    This is an experimental class for nodal analysis.  The API
    is likely to change.

    >>> from lcapy import Circuit, NodalAnalysis
    >>> cct = Circuit('''
    ... V1 1 0 {u(t)}; down
    ... R1 1 2; right=2
    ... L1 2 3; down=2
    ... W1 0 3; right
    ... W 1 5; up
    ... W 2 6; up
    ... C1 5 6; right=2
    ...''')    

    To perform nodal analysis in the Laplace domain:

    >>> na = NodalAnalysis(cct.laplace)
    
    To display the system of equations (in matrix form) that needs to
    be solved:
    
    >>> na.equations().pprint()
    
    To display the equations found by applying KCL at each node:
    
    >>> na.nodal_equations().pprint()

    """

    def __init__(self, cct):

        # In principle we could generate the system of
        # integro-differential equations in the time domain.
        if cct.kind == 'super':
            raise ValueError('Can only handle sub-circuits at present')
        
        self.cct = cct
        self.G = CircuitGraph(cct)

        self.ydict = self._make_unknowns()
        
        self.y = matrix([val for key, val in self.ydict.items() if key != '0'])
        
        self.equations_dict = self._make_equations()

        self.A, self.b = self._analyse()
        
        
    def nodes(self):

        return self.G.nodes()

    def _make_unknowns(self):

        # Determine node voltage variables.
        ydict = ExprDict()

        # cct.node_list is sorted alphabetically
        for node in self.cct.node_list:
            if node == '0':
                ydict[node] = 0
            else:
                ydict[node] = self.Vname('Vn' + node)
        return ydict

    def _make_equations(self):

        equations_dict = ExprDict()
        for node in self.nodes():
            if node == '0':
                continue

            voltage_sources = []
            for elt in self.G.connected(node):
                if elt.type == 'V':
                    voltage_sources.append(elt)

            if voltage_sources != []:
                elt = voltage_sources[0]
                n1 = self.G.node_map[elt.nodenames[0]]
                n2 = self.G.node_map[elt.nodenames[1]]

                V = Voltage(elt.voc).select(self.cct.kind)
                
                eq = equation(self.ydict[n1], self.ydict[n2] + V)

            else:
                expr = 0
                for elt in self.G.connected(node):
                    if len(elt.nodenames) < 2:
                        raise ValueError('Elt %s has too few nodes' % elt)
                    n1 = self.G.node_map[elt.nodenames[0]]
                    n2 = self.G.node_map[elt.nodenames[1]]
                    if node == n1:
                        pass
                    elif node == n2:
                        n1, n2 = n2, n1
                    else:
                        raise ValueError('Component %s does not have node %d' % (elt, node))
                    expr += (self.ydict[n1] - self.ydict[n2]) * elt.cpt.Y
                    
                eq = equation(expr, 0)

            equations_dict[node] = eq

        return equations_dict        

    def _analyse(self):

        eqns = matrix(list(self.equations_dict.values()))

        return sym.linear_eq_to_matrix(eqns.expand(), *self.y)

    def nodal_equations(self):
        """Return the equations found by applying KCL at each node.  This is a
        directory of equations keyed by the node name."""

        return self.equations_dict

    def equations(self):
        """Return the equations in matrix form."""
        
        return expr(equation(sym.MatMul(self.A, self.y), self.b))

from .expr import ExprDict, expr
from .texpr import Vt
from .voltage import Voltage
from .matrix import matrix
