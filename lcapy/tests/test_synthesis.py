from lcapy import *
from lcapy.sexpr import LaplaceDomainImpedance
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

    def test_seriesRL(self):

        n1 = R(2) + L(3)
        Z1 = n1.Z(s)
        n2 = Z1.network('seriesRL')
        Z2 = n2.Z(s)

        self.assertEqual(Z1, Z2, 'R + L')

    def test_seriesRC(self):

        n1 = R(2) + C(3)
        Z1 = n1.Z(s)
        n2 = Z1.network('seriesRC')
        Z2 = n2.Z(s)

        self.assertEqual(Z1, Z2, 'R + C')

    def test_seriesGC(self):

        n1 = G(2) + C(3)
        Z1 = n1.Z(s)
        n2 = Z1.network('seriesGC')
        Z2 = n2.Z(s)

        self.assertEqual(Z1, Z2, 'G + C')

    def test_seriesLC(self):

        n1 = L(2) + C(3)
        Z1 = n1.Z(s)
        n2 = Z1.network('seriesLC')
        Z2 = n2.Z(s)

        self.assertEqual(Z1, Z2, 'L + C')

    def test_parallelRL(self):

        n1 = R(2) | L(3)
        Z1 = n1.Z(s)
        n2 = Z1.network('parallelRL')
        Z2 = n2.Z(s)

        self.assertEqual(Z1, Z2, 'R | L')

    def test_parallelRC(self):

        n1 = R(2) | C(3)
        Z1 = n1.Z(s)
        n2 = Z1.network('parallelRC')
        Z2 = n2.Z(s)

        self.assertEqual(Z1, Z2, 'R | C')

    def test_parallelGC(self):

        n1 = G(2) | C(3)
        Z1 = n1.Z(s)
        n2 = Z1.network('parallelGC')
        Z2 = n2.Z(s)

        self.assertEqual(Z1, Z2, 'G | C')

    def test_parallelLC(self):

        n1 = L(2) | C(3)
        Z1 = n1.Z(s)
        n2 = Z1.network('parallelLC')
        Z2 = n2.Z(s)

        self.assertEqual(Z1, Z2, 'L | C')

    def test_Z1(self):

        Z1 = LaplaceDomainImpedance(
            9 * s**5 + 30 * s**3 + 24 * s) / (18 * s**4 + 36 * s**2 + 8)

        self.assertEqual(Z1, Z1.network('fosterI').Z(s), 'fosterI')
        self.assertEqual(Z1, Z1.network('fosterII').Z(s), 'fosterII')
        self.assertEqual(Z1, Z1.network('cauerI').Z(s), 'cauerI')
        self.assertEqual(Z1, Z1.network('cauerII').Z(s), 'cauerII')

    def test_Z2(self):

        # Cauer II form.
        n = ((C(2) | L(3)) + C(4)) | L(5)
        Z1 = n.Z(s)

        self.assertEqual(Z1, Z1.network('fosterI').Z(s), 'fosterI')
        self.assertEqual(Z1, Z1.network('fosterII').Z(s), 'fosterII')
        self.assertEqual(Z1, Z1.network('cauerI').Z(s), 'cauerI')
        self.assertEqual(Z1, Z1.network('cauerII').Z(s), 'cauerII')

    def test_Z3(self):

        n = (R(2) | L(3)) + (R(4) | L(5))
        Z1 = n.Z(s)

        self.assertEqual(Z1, Z1.network('fosterI').Z(s), 'fosterI')
        self.assertEqual(Z1, Z1.network('fosterII').Z(s), 'fosterII')
        self.assertEqual(Z1, Z1.network('cauerI').Z(s), 'cauerI')
        self.assertEqual(Z1, Z1.network('cauerII').Z(s), 'cauerII')
