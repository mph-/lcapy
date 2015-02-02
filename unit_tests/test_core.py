from lcapy import *
import unittest
import sympy as sym

s = sym.var('s')

class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy

    """

    def test_sExpr(self):
        """Lcapy: check sExpr

        """
        a = sExpr('(s+2)/(s-2)')
        self.assertEqual(a.N, sExpr('s+2'), "N incorrect.")
        self.assertEqual(a.D, sExpr('s-2'), "D incorrect.")

