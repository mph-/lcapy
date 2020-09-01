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

    def test_Impedance(self):

        ZC = C(3).Z

        self.assertEqual(ZC.Y, 3 * j * omega, '.Y')
        self.assertEqual(ZC.Z, -j / (3 * omega), '.Z')        
        self.assertEqual(ZC.R, 0, '.R')
        self.assertEqual(ZC.G, 0, '.G')
        self.assertEqual(ZC.X, -1 / (3 * omega), '.X')
        self.assertEqual(ZC.B, -3 * omega, '.B')

        ZR = R(3).Z

        self.assertEqual(ZR.Y, expr('1 / 3'), '.Y')
        self.assertEqual(ZR.Z, 3, '.Z')        
        self.assertEqual(ZR.R, 3, '.R')
        self.assertEqual(ZR.G, expr('1 / 3'), '.G')
        self.assertEqual(ZR.X, 0, '.X')
        self.assertEqual(ZR.B, 0, '.B')

        ZL = L(3).Z

        self.assertEqual(ZL.Y, -j / (3 * omega), '.Y')
        self.assertEqual(ZL.Z, 3 * j * omega, '.Z')        
        self.assertEqual(ZL.R, 0, '.R')
        self.assertEqual(ZL.G, 0, '.G')
        self.assertEqual(ZL.X, 3 * omega, '.X')
        self.assertEqual(ZL.B, 1 / (3 * omega), '.B')                

        Z1 = (R(2) + L(3)).Z

        self.assertEqual(Z1.Y, 1 / (3 * j * omega + 2), '.Y')
        self.assertEqual(Z1.Z, 3 * j * omega + 2, '.Z')        
        self.assertEqual(Z1.R, 2, '.R')
        self.assertEqual(Z1.G, 2 / (9 * omega**2 + 4), '.G')
        self.assertEqual(Z1.X, 3 * omega, '.X')
        self.assertEqual(Z1.B, 3 * omega / (9 * omega**2 + 4), '.B')

    def test_reciprocal(self):

        Z1 = Impedance(2)
        Y1 = Admittance(1 / 2)

        # FIXME for python 2.7
        # self.assertEqual(1 / Z1, Y1, '1 / Z')
        # self.assertEqual(1 / Y1, Z1, '1 / Y')                
