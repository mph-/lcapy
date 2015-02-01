from lcapy import Circuit, R, C
from lcapy.core import Zs
import unittest
import sympy as sym

s = sym.var('s')

class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy

    """

    def test_RL1(self):
        """Lcapy: check RL

        """
        a = Circuit()
        a.add('R1 1 2')
        a.add('C1 2 0')

        self.assertEqual(a.Z(1, 2), R('R1').Z, "Z incorrect for R1.")
        self.assertEqual(a.Z(2, 0), C('C1').Z, "Z incorrect for C1.")
        self.assertEqual(a.Z(1, 0), (R('R1') + C('C1')).Z, "Z incorrect for R1 + C1.")

