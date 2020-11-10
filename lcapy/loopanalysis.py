"""
This module performs loop analysis.  It is primarily for showing
the equations rather than evaluating them.

Copyright 2019--2020 Michael Hayes, UCECE

"""

from .circuitgraph import CircuitGraph
from .expr import expr
from .utils import equation
import sympy as sym


class LoopAnalysis(object):

    def __init__(self, cct):

        self.cct = cct
        self.G = CircuitGraph(cct)

    def loops(self):

        return self.G.loops()

    def mesh_currents(self):

        if not self.G.is_planar:
            raise ValueError('Circuit topology is not planar')

        loops = self.loops()
        Nloops = len(loops)
        
        mesh_currents = ExprList([expr('i_%d(t)' % (m + 1)) for m in range(Nloops)])
        return mesh_currents
    
    def equations(self):
        """Return loop equations as a list."""

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

            result = expr(0)

            loop1 = loop.copy()
            loop1.append(loop1[0])
            for j in range(len(loop1) - 1):

                elt = self.G.component(loop1[j], loop1[j + 1])

                if elt.is_current_source:
                    raise ValueError('TODO: handle current source in loop')
                
                reversed = elt.nodenames[0] == loop1[j]

                current = mesh_currents[m]

                # Map node names to equipotential node names.
                nodenames = [self.cct.node_map[nodename] for nodename in elt.nodenames]
                
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
                            current += mesh_currents[n]
                            break
                        elif (nodenames[1] == loop2[l] and
                            nodenames[0] == loop2[l + 1]):
                            current -= mesh_currents[n]
                            break
                            
                v = elt.cpt.v_equation(current)                
                result += v

            eq = equation(result, 0)        
            equations.append(eq)

        return equations
    

from .expr import ExprList, expr    
