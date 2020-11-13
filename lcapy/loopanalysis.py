"""
This module performs loop analysis.  It is primarily for showing
the equations rather than evaluating them.

Copyright 2019--2020 Michael Hayes, UCECE

"""

# TODO, handle current sources and dependent sources.

from .circuitgraph import CircuitGraph
from .expr import expr
from .utils import equation
import sympy as sym


class LoopAnalysis(object):
    """
    This is an experimental class for loop analysis.  The API
    is likely to change.

    >>> from lcapy import Circuit, LoopAnalysis
    >>> cct = Circuit('''
    ... V1 1 0 {u(t)}; down
    ... R1 1 2; right=2
    ... L1 2 3; down=2
    ... W1 0 3; right
    ... W 1 5; up
    ... W 2 6; up
    ... C1 5 6; right=2
    ...''')    

    To perform nodal analysis in the time domain:

    >>> la = LoopAnalysis(cct)
    
    To display the equations found by applying KVL around each mesh:
    
    >>> la.mesh_equations().pprint()

    """
    
    def __init__(self, cct):

        self.cct = cct
        self.G = CircuitGraph(cct)

        self.kind = self.cct.kind
        if self.kind == 'super':
            self.kind = 'time'        

    def loops(self):

        return self.G.loops()

    @property
    def num_loops(self):

        return len(self.G.loops)

    def mesh_currents(self):

        if not self.G.is_planar:
            raise ValueError('Circuit topology is not planar')

        loops = self.loops()
        Nloops = len(loops)
        
        mesh_currents = ExprList([Current('i_%d(t)' % (m + 1)).select(self.kind) for m in range(Nloops)])
        return mesh_currents
    
    def mesh_equations_list(self):
        """Return mesh equations as a list."""

        if not self.G.is_planar:
            raise ValueError('Circuit topology is not planar')

        if self.cct.is_mixed:
            raise ValueError('Circuit has a mixture of ac/dc/transient sources')
        
        # This is work in progress to do mesh analysis.

        loops = self.loops()
        Nloops = len(loops)

        mesh_currents = [expr('i_%d(t)' % (m + 1)) for m in range(Nloops)]
        equations = ExprList()        

        for m, loop in enumerate(loops):

            result = Voltage(0).select(self.kind)

            loop1 = loop.copy()
            loop1.append(loop1[0])
            for j in range(len(loop1) - 1):

                elt = self.G.component(loop1[j], loop1[j + 1])
                if elt is None:
                    continue

                if elt.is_current_source:
                    raise ValueError('TODO: handle current source in loop')
                
                # Map node names to equipotential node names.
                nodenames = [self.cct.node_map[nodename] for nodename in elt.nodenames]

                is_reversed = nodenames[0] == loop1[j] and nodenames[1] == loop1[j + 1]
                
                if is_reversed:
                    current = -mesh_currents[m]
                else:
                    current = mesh_currents[m]                    

                for n, loop2 in enumerate(loops):

                    if loop2 == loop:
                        continue

                    loop2 = loop2.copy()
                    loop2.append(loop2[0])

                    if nodenames[0] not in loop2 or nodenames[1] not in loop2:
                        continue

                    # Determine polarity of cpt.
                    for l in range(len(loop2) - 1):
                        if (nodenames[0] == loop2[l] and
                            nodenames[1] == loop2[l + 1]):
                            current -= mesh_currents[n]
                            break
                        elif (nodenames[1] == loop2[l] and
                            nodenames[0] == loop2[l + 1]):
                            current += mesh_currents[n]
                            break
                            
                v = elt.cpt.v_equation(current, self.kind)
                if elt.is_voltage_source and is_reversed:
                    v = -v
                result += v

            eq = equation(result, 0)        
            equations.append(eq)

        return equations

    def mesh_equations(self):
        """Return mesh equations as a dict keyed by the mesh current."""

        result = ExprDict()

        for current, equation in zip(self.mesh_currents(),
                                     self.mesh_equations_list()):
            result[current] = equation
        return result

from .expr import ExprList, ExprDict, expr    
from .current import Current
from .voltage import Voltage
