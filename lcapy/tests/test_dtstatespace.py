from lcapy import *
import numpy as np
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy discrete-time state space analysis
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

    def transfer(self):

        ss = DTStateSpace(((0, 1), (1, 1)), (1, 1), ((1, 2),), [1])

        umin = ss.minimum_energy_input(2, [5, 7], [0, 0])
        self.assertEqual(umin, [2, 3], "minimum_energy_input")

        umin = ss.minimum_energy(2, [5, 7], [0, 0])
        self.assertEqual(umin, 13, "minimum_energy")

        xfinal = ss.state_transfer([[2], [3]], xinitial=[0, 0])
        self.assertEqual(xfinal, [5, 7], "state_transfer")
