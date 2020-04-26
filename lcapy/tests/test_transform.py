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

    def test_omega_to_f(self):

        a = 3 * omega
        b = a(f)
        c = a.transform(f)
        
        self.assertEqual(b, 6 * pi * f, 'transform(f)')
        self.assertEqual(c, 6 * pi * f, 'transform(f)')

    def test_f_to_omega(self):

        a = 2 * pi * f
        b = a(omega)
        c = a.transform(omega)
        
        self.assertEqual(b, omega, 'transform(omega)')
        self.assertEqual(c, omega, 'transform(omega)')                
