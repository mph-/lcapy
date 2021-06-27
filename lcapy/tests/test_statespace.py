from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy state space analysis    
    """

    def assertEqual2(self, ans1, ans2, comment):

        ans1 = ans1.canonical()
        ans2 = ans2.canonical()

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            ans1.pprint()
            ans2.pprint()
            raise AssertionError(e)

    def test_VRC1(self):
        """Check VRC circuit

        """
        a = Circuit("""
        V1 1 0
        R1 1 2
        C1 2 0""")

        ss = a.ss
        
        self.assertEqual2(expr(ss.x[0]), expr('v_C1(t)'), "Incorrect state variable")
        self.assertEqual2(expr(ss.y[0]), expr('v_1(t)'), "Incorrect output variable1")
        self.assertEqual2(expr(ss.y[1]), expr('v_2(t)'), "Incorrect output variable2") 
        self.assertEqual2(expr(ss.A[0]), expr('-1/(R1 * C1)'), "Incorrect A matrix")
        self.assertEqual2(expr(ss.B[0]), expr('1/(R1 * C1)'), "Incorrect B matrix")
        self.assertEqual2(expr(ss.C[0]), expr(0), "Incorrect C[0] matrix element")
        self.assertEqual2(expr(ss.C[1]), expr(1), "Incorrect C[1] matrix element")
        self.assertEqual2(expr(ss.D[0]), expr(1), "Incorrect D[0] matrix element")
        self.assertEqual2(expr(ss.D[1]), expr(0), "Incorrect D[1] matrix element")
        self.assertEqual2(expr(ss.eigenvalues[0]), expr('-1/(R1 * C1)'), "Incorrect eigenvalues")
        
        
    def test_VRL1(self):
        """Check VRL circuit

        """
        a = Circuit("""
        V1 1 0
        R1 1 2
        L1 2 0""")

        ss = a.ss
        
        self.assertEqual2(expr(ss.x[0]), expr('i_L1(t)'), "Incorrect state variable")
        self.assertEqual2(expr(ss.y[0]), expr('v_1(t)'), "Incorrect output variable1")
        self.assertEqual2(expr(ss.y[1]), expr('v_2(t)'), "Incorrect output variable2")        
        self.assertEqual2(expr(ss.A[0]), expr('-R1 / L1'), "Incorrect A matrix")
        self.assertEqual2(expr(ss.B[0]), expr('1 / L1'), "Incorrect B matrix")
        self.assertEqual2(expr(ss.C[0]), expr(0), "Incorrect C[0] matrix element")
        self.assertEqual2(expr(ss.C[1]), expr('-R1'), "Incorrect C[1] matrix element")
        self.assertEqual2(expr(ss.D[0]), expr(1), "Incorrect D[0] matrix element")
        self.assertEqual2(expr(ss.D[1]), expr(1), "Incorrect D[1] matrix element")
        self.assertEqual2(expr(ss.eigenvalues[0]), expr('-R1 / L1'), "Incorrect eigenvalues")
    
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
        
