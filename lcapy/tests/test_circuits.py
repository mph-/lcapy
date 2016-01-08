from lcapy import Circuit, R, C, L, V, I, v
from lcapy.core import Zs, s
import unittest
import sympy as sym


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

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

    def test_RL1(self):
        """Lcapy: check RL network

        """
        a = Circuit()
        a.add('R1 1 2')
        a.add('C1 2 0')

        self.assertEqual2(a.Z(1, 2), R('R1').Z, "Z incorrect for R1.")
        self.assertEqual2(a.Z(2, 0), C('C1').Z, "Z incorrect for C1.")
        self.assertEqual2(
            a.Z(1, 0), (R('R1') + C('C1')).Z, "Z incorrect for R1 + C1.")

        self.assertEqual2(a.Y(1, 2), R('R1').Y, "Y incorrect for R1.")
        self.assertEqual2(a.Y(2, 0), C('C1').Y, "Y incorrect for C1.")
        self.assertEqual2(
            a.Y(1, 0), (R('R1') + C('C1')).Y, "Y incorrect for R1 + C1.")
        self.assertEqual2(a.Isc(1, 0), I(0).I, "Isc incorrect")

    def test_VRC1(self):
        """Lcapy: check VRC circuit

        """
        a = Circuit()
        a.add('V1 1 0')
        a.add('R1 1 2')
        a.add('C1 2 0')

        # Note, V1 acts as a short-circuit for the impedance/admittance
        self.assertEqual2(
            a.Z(1, 2), (R('R1') | C('C1')).Z, "Z incorrect across R1")
        self.assertEqual2(
            a.Z(2, 0), (R('R1') | C('C1')).Z, "Z incorrect across C1")
        self.assertEqual2(a.Z(1, 0), R(0).Z, "Z incorrect across V1")

        self.assertEqual2(
            a.Y(1, 2), (R('R1') | C('C1')).Y, "Y incorrect across R1")
        self.assertEqual2(
            a.Y(2, 0), (R('R1') | C('C1')).Y, "Y incorrect across C1")
        # This has a non-invertible A matrix.
        # self.assertEqual2(a.Y(1, 0), R(0).Y, "Y incorrect across V1")

        self.assertEqual2(a.Voc(1, 0), V('V1').V, "Voc incorrect across V1")


    def test_VRL1(self):
        """Lcapy: check VRL circuit

        """
        a = Circuit()
        a.add('V1 1 0')
        a.add('R1 1 2')
        a.add('L1 2 0')

        # Note, V1 acts as a short-circuit for the impedance/admittance
        self.assertEqual2(
            a.Z(1, 2), (R('R1') | L('L1')).Z, "Z incorrect across R1")
        self.assertEqual2(
            a.Z(2, 0), (R('R1') | L('L1')).Z, "Z incorrect across L1")
        self.assertEqual2(a.Z(1, 0), R(0).Z, "Z incorrect across V1")

        self.assertEqual2(
            a.Y(1, 2), (R('R1') | L('L1')).Y, "Y incorrect across R1")
        self.assertEqual2(
            a.Y(2, 0), (R('R1') | L('L1')).Y, "Y incorrect across L1")
        # This has a non-invertible A matrix.
        # self.assertEqual2(a.Y(1, 0), R(0).Y, "Y incorrect across V1")

        self.assertEqual2(a.Voc(1, 0), V('V1').V, "Voc incorrect across V1")


    def test_IR1(self):
        """Lcapy: check IR circuit

        """
        a = Circuit()
        a.add('I1 1 0 2')
        a.add('R1 1 0 1')

        self.assertEqual2(a.R1.V, V(2).V, "Incorrect voltage")

        self.assertEqual2(a[1].V, V(2).V, "Incorrect node voltage")

    def test_VCVS1(self):
        """Lcapy: check VCVS

        """
        a = Circuit()
        a.add('V1 1 0 2')
        a.add('R1 1 0 1')
        a.add('E1 2 0 1 0 3')
        a.add('R2 2 0 1')

        self.assertEqual2(a.R2.V, V(6).V, "Incorrect voltage")

    def test_VCCS1(self):
        """Lcapy: check VCCS

        """
        a = Circuit()
        a.add('V1 1 0 2')
        a.add('R1 1 0 1')
        a.add('G1 2 0 1 0 3')
        a.add('R2 2 0 1')

        self.assertEqual2(a.R2.V, V(6).V, "Incorrect voltage")

    def test_CCCS1(self):
        """Lcapy: check CCCS

        """
        a = Circuit()
        a.add('V1 1 0 10')
        a.add('R1 1 2 2')
        a.add('V2 2 0 0')
        a.add('F1 3 0 V2 2')
        a.add('R2 3 0 1')

        self.assertEqual2(a.R2.V, V(10).V, "Incorrect voltage")


    def test_CCVS1(self):
        """Lcapy: check CCVS

        """
        a = Circuit()
        a.add('V1 1 0 10')
        a.add('R1 1 2 2')
        a.add('V2 2 0 0')
        a.add('H1 3 0 V2 2')
        a.add('R2 3 0 1')

        self.assertEqual2(a.R2.V, V(10).V, "Incorrect voltage")


    def test_V1(self):
        """Lcapy: test V1"""

        a = Circuit()
        a.add('V1 1 0 10') 

        self.assertEqual2(a.V1.V, 10 / s, "Incorrect voltage")
