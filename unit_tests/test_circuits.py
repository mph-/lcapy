from lcapy import Circuit, R, C
from lcapy.core import Zs
import unittest
import sympy as sym

s = sym.var('s')

class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy

    """

    def test_RL1(self):
        """Lcapy: check RL network

        """
        a = Circuit()
        a.add('R1 1 2')
        a.add('C1 2 0')

        self.assertEqual(a.Z(1, 2), R('R1').Z, "Z incorrect for R1.")
        self.assertEqual(a.Z(2, 0), C('C1').Z, "Z incorrect for C1.")
        self.assertEqual(a.Z(1, 0), (R('R1') + C('C1')).Z, "Z incorrect for R1 + C1.")

        self.assertEqual(a.Y(1, 2), R('R1').Y, "Y incorrect for R1.")
        self.assertEqual(a.Y(2, 0), C('C1').Y, "Y incorrect for C1.")
        self.assertEqual(a.Y(1, 0).canonical(), (R('R1') + C('C1')).Y.canonical(), "Y incorrect for R1 + C1.")


    def test_VRL1(self):
        """Lcapy: check RL circuit

        """
        a = Circuit()
        a.add('V1 1 0')
        a.add('R1 1 2')
        a.add('C1 2 0')

        # Note, V1 acts as a short-circuit for the impedance/admittanc
        self.assertEqual(a.Z(1, 2).canonical(), (R('R1') | C('C1')).Z.canonical(), "Z incorrect across R1")
        self.assertEqual(a.Z(2, 0).canonical(), (R('R1') | C('C1')).Z.canonical(), "Z incorrect across C1")
        self.assertEqual(a.Z(1, 0).canonical(), R(0).Z.canonical(), "Z incorrect across V1")

        self.assertEqual(a.Y(1, 2).canonical(), (R('R1') | C('C1')).Y.canonical(), "Y incorrect across R1")
        self.assertEqual(a.Y(2, 0).canonical(), (R('R1') | C('C1')).Y.canonical(), "Y incorrect across C1")
        # This has a non-invertible A matrix.
        # self.assertEqual(a.Y(1, 0).canonical(), R(0).Y.canonical(), "Y incorrect across V1")
