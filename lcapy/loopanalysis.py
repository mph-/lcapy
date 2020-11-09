"""
This module performs loop analysis.  It is primarily for showing
the equations rather than evaluating them.

Copyright 2019--2020 Michael Hayes, UCECE

"""

from .circuitgraph import CircuitGraph


class LoopAnalysis(object):

    def __init__(self, cct):

        self.cct = cct
        self.G = CircuitGraph(cct)

    def loops(self):

        return self.G.loops()

    def analyse(self):

        loops = self.loops()
        Nloops = len(loops)

        mesh_currents = ['i_%d' % (m + 1) for m in range(Nloops)]

        for m, loop in enumerate(loops):

            print(loop)
            loop.append(loop[0])
            for j in range(len(loop) - 1):

                elt = self.G.component(loop[j], loop[j + 1])
                reversed = elt.nodenames[0] == loop[j]

                print(elt, reversed, elt.is_reactance)
                
        


    
