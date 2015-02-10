from lcapy import *
import unittest

class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy

    """

    def assertEqual2(self, ans1, ans2, comment):

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError, e:
            ans1.pprint()
            ans2.pprint()
            raise AssertionError(e)


    def test_sExpr1(self):
        """Lcapy: check sExpr

        """
        a = sExpr('(s+2)/(s-2)')
        self.assertEqual2(a.N, sExpr('s+2'), "N incorrect.")
        self.assertEqual2(a.D, sExpr('s-2'), "D incorrect.")

        self.assertEqual2(a.poles().keys(), [2], "poles incorrect.")
        self.assertEqual2(a.zeros().keys(), [-2], "zeros incorrect.")

        self.assertEqual2(a.partfrac(), 1 + 4 / (s - 2), "partfrac incorrect.")


    def test_sExpr2(self):
        """Lcapy: check sExpr 2

        """
        a = (s + 2) * (s + 3) / (s - 2)
        self.assertEqual2(a.N, (s + 2) * (s + 3), "N incorrect.")
        self.assertEqual2(a.D, s - 2, "D incorrect.")

        self.assertEqual2(a.poles().keys(), [2], "poles incorrect.")
        self.assertEqual2(a.zeros().keys(), [-2, -3], "zeros incorrect.")

        self.assertEqual2(a.partfrac(), s + 7 + 20 / (s - 2), "partfrac incorrect.")
        self.assertEqual2(a.mixedfrac(), s + 7 + 20 / (s - 2), "mixedfrac incorrect.")
        self.assertEqual2(a.general(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")
        self.assertEqual2(a.canonical(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")


    def test_sExpr3(self):
        """Lcapy: check sExpr 3

        """
        a = (s**2 + 5 * s + 6) / (s - 2)
        self.assertEqual2(a.N, s ** 2 + 5 * s + 6, "N incorrect.")
        self.assertEqual2(a.D, s - 2, "D incorrect.")

        self.assertEqual2(a.poles().keys(), [2], "poles incorrect.")
        self.assertEqual2(a.zeros().keys(), [-2, -3], "zeros incorrect.")

        self.assertEqual2(a.partfrac(), s + 7 + 20 / (s - 2), "partfrac incorrect.")
        self.assertEqual2(a.mixedfrac(), s + 7 + 20 / (s - 2), "mixedfrac incorrect.")
        self.assertEqual2(a.general(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")
        self.assertEqual2(a.canonical(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")
        self.assertEqual2(a.ZPK(), (s + 2) * (s + 3) / (s - 2), "ZPK incorrect.")


