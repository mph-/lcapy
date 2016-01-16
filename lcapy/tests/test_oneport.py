from lcapy import *
import unittest
import sympy as sym

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

    def test_R(self):
        """Lcapy: check R

        """
        a = R(10)
        self.assertEqual2(a.Z, 10, "Z incorrect.")

    def test_L(self):
        """Lcapy: check L

        """
        a = L(10, 5)
        self.assertEqual2(a.Z, 10 * s, "Z incorrect.")
        self.assertEqual2(a.V, -10 * 5, "V incorrect.")

    def test_C(self):
        """Lcapy: check C

        """
        a = C(10, 5)

        self.assertEqual2(a.Z, 1 / (10 * s), "Z incorrect.")
        self.assertEqual2(a.V, 5 / s, "V incorrect.")

    def test_R_series_R(self):
        """Lcapy: check R + R

        """
        a = R(10) + R(5)
        b = a.simplify()
        self.assertEqual2(b.R, 15, "R incorrect.")
        self.assertEqual2(b.Z, 15, "Z incorrect.")
        self.assertEqual2(type(b), R, "type incorrect.")

    def test_L_series_L(self):
        """Lcapy: check L + L

        """
        a = L(10) + L(5)
        b = a.simplify()
        self.assertEqual2(b.L, 15, "L incorrect.")
        self.assertEqual2(b.Z, 15 * s, "Z incorrect.")
        self.assertEqual2(type(b), L, "type incorrect.")

    def test_C_series_C(self):
        """Lcapy: check C + C

        """
        a = C(10) + C(15)
        b = a.simplify()
        self.assertEqual2(b.C, 6, "C incorrect.")
        self.assertEqual2(b.Z, 1 / (6 * s), "Z incorrect.")
        self.assertEqual2(type(b), C, "type incorrect.")

    def test_V_series_V(self):
        """Lcapy: check Vdc + Vdc

        """
        a = Vdc(10) + Vdc(5)
        b = a.simplify()
        self.assertEqual2(b.v, Vdc(15).v, "V incorrect.")
        self.assertEqual2(b.V, 15 / s, "V incorrect.")
        self.assertEqual2(type(b), Vdc, "type incorrect.")

    def test_R_series_L(self):
        """Lcapy: check R + L

        """
        a = R(10) + L(5)
        self.assertEqual2(a.Z, 5 * s + 10, "Z incorrect.")

    def test_R_parallel_R(self):
        """Lcapy: check R | R

        """
        a = R(10) | R(15)
        b = a.simplify()
        self.assertEqual2(b.R, 6, "R incorrect.")
        self.assertEqual2(b.Z, 6, "Z incorrect.")
        self.assertEqual2(type(b), R, "type incorrect.")

    def test_L_parallel_L(self):
        """Lcapy: check L | L

        """
        a = L(10) | L(15)
        b = a.simplify()
        self.assertEqual2(b.L, 6, "L incorrect.")
        self.assertEqual2(b.Z, 6 * s, "Z incorrect.")
        self.assertEqual2(type(b), L, "type incorrect.")

    def test_C_parallel_C(self):
        """Lcapy: check C | C

        """
        a = C(10) | C(15)
        b = a.simplify()
        self.assertEqual2(b.C, 25, "C incorrect.")
        self.assertEqual2(b.Z, 1 / (25 * s), "Z incorrect.")
        self.assertEqual2(type(b), C, "type incorrect.")

    def test_I_parallel_I(self):
        """Lcapy: check I | I

        """
        a = Idc(10) | Idc(5)
        b = a.simplify()
        self.assertEqual2(b.i, Idc(15).i, "I incorrect.")
        self.assertEqual2(b.I, 15 / s, "I incorrect.")
        self.assertEqual2(type(b), Idc, "type incorrect.")

    def test_load(self):
        """Lcapy: check load

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.load(R(30))

        self.assertEqual2(type(b), Thevenin, "type incorrect.")
        self.assertEqual2(b.Z, 7.5, "Shunt loaded R incorrect Z.")

    def test_open_circuit(self):
        """Lcapy: check open_circuit

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.open_circuit()

        self.assertEqual2(type(b), Thevenin, "type incorrect.")
        self.assertEqual2(b.Z, 10, "incorrect Z.")
        self.assertEqual2(b.V, 5 / s, "incorrect V.")

    def test_short_circuit(self):
        """Lcapy: check short_circuit

        """

        a = Series(R(10) + Vdc(5))
        b = a.short_circuit()

        self.assertEqual2(type(b), Norton, "type incorrect.")
        self.assertEqual2(b.Z, 10, "incorrect Z.")
        self.assertEqual2(b.V, 5 / s, "incorrect V.")
