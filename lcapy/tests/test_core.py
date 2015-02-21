from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def assertEqual2(self, ans1, ans2, comment):

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            ans1.pprint()
            ans2.pprint()
            raise AssertionError(e)

    def test_sExpr1(self):
        """Lcapy: check sExpr1

        """
        a = sExpr('(s+2)/(s-2)')
        self.assertEqual2(a.N, sExpr('s+2'), "N incorrect.")
        self.assertEqual2(a.D, sExpr('s-2'), "D incorrect.")

        self.assertEqual2(a.poles().keys(), [2], "poles incorrect.")
        self.assertEqual2(a.zeros().keys(), [-2], "zeros incorrect.")

        self.assertEqual2(a.partfrac(), 1 + 4 / (s - 2), "partfrac incorrect.")

        self.assertEqual(a.evaluate(1), -3.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate(-2), 0.0, "scalar evaluate incorrect.")

        self.assertEqual(a.inverse_laplace(), 4 * exp(2 * t) * H(t) + DiracDelta(t), "inverse Laplace incorrect.")

        aw = a(j * omega)
        self.assertEqual2(aw.real, (omega**2 - 4) / (omega**2 + 4), "real part incorrect.")
        self.assertEqual2(aw.imag, -4 * omega / (omega**2 + 4), "imag part incorrect.")


    def test_sExpr2(self):
        """Lcapy: check sExpr2

        """
        a = (s + 2) * (s + 3) / (s - 2)
        self.assertEqual2(a.N, (s + 2) * (s + 3), "N incorrect.")
        self.assertEqual2(a.D, s - 2, "D incorrect.")

        self.assertEqual2(a.poles().keys(), [2], "poles incorrect.")
        self.assertEqual2(a.zeros().keys(), [-2, -3], "zeros incorrect.")

        self.assertEqual2(
            a.partfrac(), s + 7 + 20 / (s - 2), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), s + 7 + 20 / (s - 2), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")
        self.assertEqual2(
            a.canonical(), (s**2 + 5 * s + 6) / (s - 2), "canonical incorrect.")

        self.assertEqual(a.inverse_laplace(), 20 * exp(2 * t) * H(t) + 7 * DiracDelta(t) + DiracDelta(t, 1), "inverse Laplace incorrect.")

    def test_sExpr3(self):
        """Lcapy: check sExpr3

        """
        a = (s**2 + 5 * s + 6) / (s - 2)
        self.assertEqual2(a.N, s ** 2 + 5 * s + 6, "N incorrect.")
        self.assertEqual2(a.D, s - 2, "D incorrect.")

        self.assertEqual2(a.poles().keys(), [2], "poles incorrect.")
        self.assertEqual2(a.zeros().keys(), [-2, -3], "zeros incorrect.")

        self.assertEqual2(
            a.partfrac(), s + 7 + 20 / (s - 2), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), s + 7 + 20 / (s - 2), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")
        self.assertEqual2(
            a.canonical(), (s**2 + 5 * s + 6) / (s - 2), "canonical incorrect.")
        self.assertEqual2(
            a.ZPK(), (s + 2) * (s + 3) / (s - 2), "ZPK incorrect.")

        self.assertEqual(a.inverse_laplace(), 20 * exp(2 * t) * H(t) + 7 * DiracDelta(t) + DiracDelta(t, 1), "inverse Laplace incorrect.")


    def test_sExpr4(self):
        """Lcapy: check sExpr4

        """
        a = 1 / ((s - j) * (s + j))
        self.assertEqual2(a.N, 1, "N incorrect.")
        self.assertEqual2(a.D, (s - j) * (s + j), "D incorrect.")

        self.assertEqual2(a.poles().keys(), [-j, j], "poles incorrect.")
        self.assertEqual2(a.zeros().keys(), [], "zeros incorrect.")

        self.assertEqual2(
            a.partfrac(), 0.5j / (s + j) - 0.5j / (s - j), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), 1 / (s**2 + 1), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), 1 / (s**2 + 1), "general incorrect.")
        self.assertEqual2(
            a.canonical(), 1 / (s**2 + 1), "canonical incorrect.")
        self.assertEqual2(
            a.ZPK(), 1 / ((s - j) * (s + j)), "ZPK incorrect.")

        self.assertEqual(a.inverse_laplace(), -sin(t) * H(t), "inverse Laplace incorrect.")


    def test_sExpr5(self):
        """Lcapy: check sExpr5

        """
        a = 1 / ((s - 1j) * (s + 1j))
        self.assertEqual2(a.N, 1, "N incorrect.")
        self.assertEqual2(a.D, (s - j) * (s + j), "D incorrect.")

        self.assertEqual2(a.poles().keys(), [-j, j], "poles incorrect.")
        self.assertEqual2(a.zeros().keys(), [], "zeros incorrect.")

        self.assertEqual2(
            a.partfrac(), 0.5j / (s + j) - 0.5j / (s - j), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), 1 / (s**2 + 1), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), 1 / (s**2 + 1), "general incorrect.")
        self.assertEqual2(
            a.canonical(), 1 / (s**2 + 1), "canonical incorrect.")
        self.assertEqual2(
            a.ZPK(), 1 / ((s - j) * (s + j)), "ZPK incorrect.")

        self.assertEqual(a.inverse_laplace(), -sin(t) * H(t), "inverse Laplace incorrect.")



    def test_tExpr1(self):
        """Lcapy: check tExpr1

        """
        a = t**2

        self.assertEqual(a.evaluate(2), 4.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate((2, 3))[1], 9.0, "vector evaluate incorrect.")


    def test_step(self):
        """Lcapy: check step

        """
        a = u(t)

        self.assertEqual(a.evaluate(2), 1.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate((-2, 2))[1], 1.0, "vector evaluate incorrect.")
        self.assertEqual(a.evaluate((-2, 2))[0], 0.0, "vector evaluate incorrect.")
        self.assertEqual(a.laplace(), 1 / s, "Laplace transform incorrect.")


    def test_delta(self):
        """Lcapy: check delta

        """
        a = delta(t)

        self.assertEqual(a.evaluate(2), 0.0, "scalar evaluate incorrect.")
        # Well, I reckon that the Laplace transfrom should integrate
        # from 0- to ensure symmetry, but at the moment 0.5 is 
        # the correct answer.
        self.assertEqual(a.laplace(), 0.5, "Laplace transform incorrect.")

    
    def test_jomega(self):
        """Lcapy: check jomega

        """

        a = sExpr('s+2')
        b = a(j * omega)

        self.assertEqual2(b, j * omega + 2, "Substitution failed.")
        self.assertEqual2(a.jomega(), j * omega + 2, "jomega failed.")


    def test_subs1(self):
        """Lcapy: check subs

        """

        a = s(omega)

        self.assertEqual(a.expr.is_real, True, "Lost is_real.")
        self.assertEqual2(a, omega, "Substitution fail.")
