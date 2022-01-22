from lcapy import *
from lcapy.sym import tausym
from lcapy.circuitgraph import CircuitGraph
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_circuitgraph(self):
        """Test circuitgraph"""

        a = Circuit('lcapy/tests/test2.sch')
        cg = CircuitGraph.from_circuit(a)

        self.assertEqual(cg.in_parallel('L'), set({'L', 'R2'}), "in_parallel")
        self.assertEqual(cg.in_series('C'), set({'C', 'R', 'V'}), "in_series")

        self.assertEqual(cg.connected('2'), set({'V', 'R'}), "connected")

        loops = cg.loops()
        loop1 = loops[0]
        loop2 = loops[1]
        if len(loop1) < len(loop2):
            loop1, loop2 = loop2, loop1

        self.assertEqual(loop1, ['0', '2', '1', '3'], "loop1")
        self.assertEqual(loop2, ['*0', '0', '3'], "loop2")

        self.assertEqual(cg.component('1', '2').name, 'R', "connected")

        self.assertEqual(cg.is_connected, True, "is_connected")
        self.assertEqual(cg.is_planar, True, "is_planar")

        self.assertEqual(cg.num_parts, 1, "num_parts")
        self.assertEqual(cg.num_nodes, 5, "num_nodes")
        self.assertEqual(cg.num_branches, 6, "num_branches")
        self.assertEqual(cg.rank, 4, "rank")
        self.assertEqual(cg.nullity, 2, "nullity")
        self.assertEqual(cg.node_connectivity, 2, "node_connectivity")

        ct = cg.tree()

        self.assertEqual(ct.num_parts, 1, "num_parts")
        self.assertEqual(ct.num_nodes, 5, "num_nodes")
        self.assertEqual(ct.num_branches, 4, "num_branches")
        self.assertEqual(ct.rank, 4, "rank")
        self.assertEqual(ct.nullity, 0, "nullity")
        self.assertEqual(ct.node_connectivity, 1, "node_connectivity")

        cl = cg.links()

        self.assertEqual(cl.num_parts, 1, "num_parts")
        self.assertEqual(cl.num_nodes, 3, "num_nodes")
        self.assertEqual(cl.num_branches, 2, "num_branches")
        self.assertEqual(cl.rank, 2, "rank")
        self.assertEqual(cl.nullity, 0, "nullity")
        self.assertEqual(cl.node_connectivity, 1, "node_connectivity")
