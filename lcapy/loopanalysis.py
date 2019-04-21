"""
This module performs loop analysis.  It is primarily for showing
the equations rather than evaluating them.

Copyright 2019 Michael Hayes, UCECE

"""

from .circuitgraph import CircuitGraph


class LoopAnalysis(object):

    def __init__(self, cct):

        self.cct = cct
        self.G = CircuitGraph(cct)

    def loops(self):

        return self.G.loops()


    def analyse(self):

        print(self.loops())
        
