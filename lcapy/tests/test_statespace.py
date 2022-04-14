from lcapy import *
import numpy as np
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

        self.assertEqual2(expr(ss.x[0]), expr(
            'v_C1(t)'), "Incorrect state variable")
        self.assertEqual2(expr(ss.y[0]), expr(
            'v_1(t)'), "Incorrect output variable1")
        self.assertEqual2(expr(ss.y[1]), expr(
            'v_2(t)'), "Incorrect output variable2")
        self.assertEqual2(expr(ss.A[0]), expr(
            '-1/(R1 * C1)'), "Incorrect A matrix")
        self.assertEqual2(expr(ss.B[0]), expr(
            '1/(R1 * C1)'), "Incorrect B matrix")
        self.assertEqual2(expr(ss.C[0]), expr(
            0), "Incorrect C[0] matrix element")
        self.assertEqual2(expr(ss.C[1]), expr(
            1), "Incorrect C[1] matrix element")
        self.assertEqual2(expr(ss.D[0]), expr(
            1), "Incorrect D[0] matrix element")
        self.assertEqual2(expr(ss.D[1]), expr(
            0), "Incorrect D[1] matrix element")
        self.assertEqual2(expr(ss.eigenvalues[0]), expr(
            '-1/(R1 * C1)'), "Incorrect eigenvalues")

    def test_VRL1(self):
        """Check VRL circuit

        """
        a = Circuit("""
        V1 1 0
        R1 1 2
        L1 2 0""")

        ss = a.ss

        self.assertEqual2(expr(ss.x[0]), expr(
            'i_L1(t)'), "Incorrect state variable")
        self.assertEqual2(expr(ss.y[0]), expr(
            'v_1(t)'), "Incorrect output variable1")
        self.assertEqual2(expr(ss.y[1]), expr(
            'v_2(t)'), "Incorrect output variable2")
        self.assertEqual2(expr(ss.A[0]), expr(
            '-R1 / L1'), "Incorrect A matrix")
        self.assertEqual2(expr(ss.B[0]), expr('1 / L1'), "Incorrect B matrix")
        self.assertEqual2(expr(ss.C[0]), expr(
            0), "Incorrect C[0] matrix element")
        self.assertEqual2(expr(ss.C[1]), expr(
            '-R1'), "Incorrect C[1] matrix element")
        self.assertEqual2(expr(ss.D[0]), expr(
            1), "Incorrect D[0] matrix element")
        self.assertEqual2(expr(ss.D[1]), expr(
            1), "Incorrect D[1] matrix element")
        self.assertEqual2(expr(ss.eigenvalues[0]), expr(
            '-R1 / L1'), "Incorrect eigenvalues")

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

    def test_transfer(self):
        Z = (s**2 + 3) / (s**3 + 2 * s + 10)

        ss = Z.state_space()

        self.assertEqual(ss.Nx, 3, "Nx")
        self.assertEqual(ss.Ny, 1, "Ny")
        self.assertEqual(ss.Nu, 1, "Nu")
        self.assertEqual(ss.is_stable, False, "is_stable")
        self.assertEqual(ss.is_controllable, True, "is_controllable")
        self.assertEqual(ss.is_observable, True, "is_observable")
        self.assertEqual(ss.is_symbolic, False, "is_symbolic")
        self.assertEqual(ss.G[0], Z, "G")

        sso = Z.state_space('OCF')
        self.assertEqual(sso.G[0], Z, "G")

        ssd = Z.state_space('DCF')
        #self.assertEqual(ssd.G[0], Z, "G")

    def test_balanced(self):

        A = [[1, -2], [3, -4]]
        B = [5, 7]
        C = [[6, 8]]
        D = [9]

        ss = StateSpace(A, B, C, D)

        ssb = ss.balance()

        h = ss.hankel_singular_values

        H = np.diag(h.numpy.squeeze())

        self.assertEqual(np.allclose(
            ssb.Wo.numpy, ssb.Wc.numpy), True, "Wo==Wc")
        self.assertEqual(np.allclose(ssb.Wo.numpy, H),
                         True, "Hankel singular values")
        self.assertEqual(ss.eigenvalues, [-1, -2], "eigen values")
