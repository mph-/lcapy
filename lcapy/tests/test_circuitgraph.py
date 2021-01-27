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

        self.assertEqual(cg.connected('1'), set({'V', 'R'}), "connected")
        
        loops = cg.loops()
        loop1 = loops[0]
        loop2 = loops[1]
        if len(loop1) < len(loop2):
            loop1, loop2 = loop2, loop1
        
        self.assertEqual(loop1, ['0', '1', '2', '3'], "loop1")
        self.assertEqual(loop2, ['*0', '0', '3'], "loop2")        

        self.assertEqual(cg.component('1', '2').name, 'R', "connected")
        
        

        
