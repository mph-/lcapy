from lcapy import *
from lcapy.sym import tausym
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_circuitgraph(self):

        a = Circuit('lcapy/tests/test2.sch')
        cg = CircuitGraph(a)

        self.assertEqual(cg.in_parallel('L'), set({'L', 'R2'}), "in_parallel")
        self.assertEqual(cg.in_series('C'), set({'C', 'R', 'V'}), "in_series")

        loops = cg.loops()
        self.assertEqual(loops, [['0', '1', '2', '3'], ['*0', '0', '3']], "loops")

        self.assertEqual(cg.component('1', '2').name, 'R', "connected")
        
        

        
