from lcapy import *
import unittest
import sympy as sym
from msignal.mrf import MRF

s = sym.var('s')


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_R(self):
        """Lcapy: check R

        """
        a = R(10)
        self.assertEqual(a.Z, 10, "Z incorrect.")

    def test_L(self):
        """Lcapy: check L

        """
        a = L(10, 5)
        self.assertEqual(a.Z, 10 * s, "Z incorrect.")
        self.assertEqual(a.V, -10 * 5, "V incorrect.")

    def test_C(self):
        """Lcapy: check C

        """
        a = C(10, 5)

        self.assertEqual(a.Z, 1 / (10 * s), "Z incorrect.")
        self.assertEqual(a.V, 5 / s, "V incorrect.")

    def test_R_series_R(self):
        """Lcapy: check R + R

        """
        a = R(10) + R(5)
        b = a.simplify()
        self.assertEqual(b.R, 15, "R incorrect.")
        self.assertEqual(b.Z, 15, "Z incorrect.")
        self.assertEqual(type(b), R, "type incorrect.")

    def test_L_series_L(self):
        """Lcapy: check L + L

        """
        a = L(10) + L(5)
        b = a.simplify()
        self.assertEqual(b.L, 15, "L incorrect.")
        self.assertEqual(b.Z, 15 * s, "Z incorrect.")
        self.assertEqual(type(b), L, "type incorrect.")

    def test_C_series_C(self):
        """Lcapy: check C + C

        """
        a = C(10) + C(15)
        b = a.simplify()
        self.assertEqual(b.C, 6, "C incorrect.")
        self.assertEqual(b.Z, 1 / (6 * s), "Z incorrect.")
        self.assertEqual(type(b), C, "type incorrect.")

    def test_V_series_V(self):
        """Lcapy: check Vdc + Vdc

        """
        a = Vdc(10) + Vdc(5)
        b = a.simplify()
        self.assertEqual(b.v, 15, "V incorrect.")
        self.assertEqual(b.V, 15 / s, "V incorrect.")
        self.assertEqual(type(b), Vdc, "type incorrect.")

    def test_R_series_L(self):
        """Lcapy: check R + L

        """
        a = R(10) + L(5)
        self.assertEqual(a.Z, 5 * s + 10, "Z incorrect.")

    def test_R_parallel_R(self):
        """Lcapy: check R | R

        """
        a = R(10) | R(15)
        b = a.simplify()
        self.assertEqual(b.R, 6, "R incorrect.")
        self.assertEqual(b.Z, 6, "Z incorrect.")
        self.assertEqual(type(b), R, "type incorrect.")

    def test_L_parallel_L(self):
        """Lcapy: check L | L

        """
        a = L(10) | L(15)
        b = a.simplify()
        self.assertEqual(b.L, 6, "L incorrect.")
        self.assertEqual(b.Z, 6 * s, "Z incorrect.")
        self.assertEqual(type(b), L, "type incorrect.")

    def test_C_parallel_C(self):
        """Lcapy: check C | C

        """
        a = C(10) | C(15)
        b = a.simplify()
        self.assertEqual(b.C, 25, "C incorrect.")
        self.assertEqual(b.Z, 1 / (25 * s), "Z incorrect.")
        self.assertEqual(type(b), C, "type incorrect.")

    def test_I_parallel_I(self):
        """Lcapy: check I | I

        """
        a = Idc(10) | Idc(5)
        b = a.simplify()
        self.assertEqual(b.i, 15, "I incorrect.")
        self.assertEqual(b.I, 15 / s, "I incorrect.")
        self.assertEqual(type(b), Idc, "type incorrect.")

    def test_load(self):
        """Lcapy: check load

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.load(R(30))

        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z, 7.5, "Shunt loaded R incorrect Z.")

    def test_open_circuit(self):
        """Lcapy: check open_circuit

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.open_circuit()

        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z, 10, "incorrect Z.")
        self.assertEqual(b.V, 5 / s, "incorrect V.")

    def test_short_circuit(self):
        """Lcapy: check short_circuit

        """

        a = Series(R(10) + Vdc(5))
        b = a.short_circuit()

        self.assertEqual(type(b), Norton, "type incorrect.")
        self.assertEqual(b.Z, 10, "incorrect Z.")
        self.assertEqual(b.V, 5 / s, "incorrect V.")
