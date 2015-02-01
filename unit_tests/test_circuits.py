from lcapy import Circuit, R, C
from lcapy.core import Zs
import unittest
import sympy as sym

s = sym.var('s')

class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy

    """

    def assertEqual2(self, ans1, ans2, comment):

        self.assertEqual(ans1.canonical(), ans2.canonical(), comment)


    def test_RL1(self):
        """Lcapy: check RL network

        """
        a = Circuit()
        a.add('R1 1 2')
        a.add('C1 2 0')

        self.assertEqual2(a.Z(1, 2), R('R1').Z, "Z incorrect for R1.")
        self.assertEqual2(a.Z(2, 0), C('C1').Z, "Z incorrect for C1.")
        self.assertEqual2(a.Z(1, 0), (R('R1') + C('C1')).Z, "Z incorrect for R1 + C1.")

        self.assertEqual2(a.Y(1, 2), R('R1').Y, "Y incorrect for R1.")
        self.assertEqual2(a.Y(2, 0), C('C1').Y, "Y incorrect for C1.")
        self.assertEqual2(a.Y(1, 0), (R('R1') + C('C1')).Y, "Y incorrect for R1 + C1.")


    def test_VRL1(self):
        """Lcapy: check RL circuit

        """
        a = Circuit()
        a.add('V1 1 0')
        a.add('R1 1 2')
        a.add('C1 2 0')

        # Note, V1 acts as a short-circuit for the impedance/admittanc
        self.assertEqual2(a.Z(1, 2), (R('R1') | C('C1')).Z, "Z incorrect across R1")
        self.assertEqual2(a.Z(2, 0), (R('R1') | C('C1')).Z, "Z incorrect across C1")
        self.assertEqual2(a.Z(1, 0), R(0).Z, "Z incorrect across V1")

        self.assertEqual2(a.Y(1, 2), (R('R1') | C('C1')).Y, "Y incorrect across R1")
        self.assertEqual2(a.Y(2, 0), (R('R1') | C('C1')).Y, "Y incorrect across C1")
        # This has a non-invertible A matrix.
        # self.assertEqual2(a.Y(1, 0), R(0).Y, "Y incorrect across V1")
