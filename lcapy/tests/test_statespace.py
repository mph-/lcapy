from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_RLC(self):

        a = Circuit("""
        V 1 0 {v(t)}; down
        R1 1 2; right
        L 2 3; right=1.5, i={i_L}
        R2 3 0_3; down=1.5, i={i_{R2}}, v={v_{R2}}
        W 0 0_3; right
        W 3 3_a; right
        C 3_a 0_4; down, i={i_C}, v={v_C}
        W 0_3 0_4; right""")
        ss = a.ss

        self.assertEqual(ss.x[0], expr('i_L(t)'), "x[0]")        
        self.assertEqual(ss.x[1], expr('v_C(t)'), "x[1]")

        self.assertEqual(ss.x0[0], 0, "x0[0]")        
        self.assertEqual(ss.x0[1], 0, "x0[1]")

        self.assertEqual(ss.y[0], expr('v_1(t)'), "y[0]")        

        self.assertEqual(ss.u[0], expr('v(t)'), "u[0]")

        self.assertEqual(ss.A[0, 0], expr('-R1/L'), "A[0, 0]")                
        
