from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_loop(self):

        a = Circuit("""
        V 1 0 {v(t)}
        R1 1 2
        R2 2 3
        R3 3 0_3
        W 0 0_3
        W 3 3_a
        R4 3_a 0_4
        W 0_3 0_4""")

        la = LoopAnalysis(a)

        la_eqs = la.mesh_equations()
        key = list(la_eqs.keys())[0]

        eq = la_eqs[key]
        # The loop ordering is random...
        if not eq.has(expr('R3')):
            key = list(la_eqs.keys())[1]            
            eq = la_eqs[key]            

        self.assertEqual(eq.lhs, voltage('R3*(I_1(t) - I_2(t)) + R4*(I_1(t) - I_2(t))'), 'mesh_equations()[1].lhs')
        self.assertEqual(eq.rhs, voltage(0), 'mesh_equations()[1].rhs')
        
        
