from lcapy import *
from lcapy.discretetime import *
import unittest
import sympy as sym


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_ztransform(self):

        self.assertEqual(nexpr(1).ZT(), z / (z - 1), "1")
        self.assertEqual(nexpr(3).ZT(), 3 * z / (z - 1), "3")
        self.assertEqual(Heaviside(n).ZT(), z / (z - 1), "Heaviside(n)")
        self.assertEqual(3 * Heaviside(n).ZT(), 3 * z / (z - 1), "3 * Heaviside(n)")
        self.assertEqual(Heaviside(n + 1).ZT(), z **2 / (z - 1), "Heaviside(n + 1)")
        self.assertEqual(Heaviside(n - 1).ZT(), 1 / (z - 1), "Heaviside(n - 1)")        
        self.assertEqual(unitimpulse(n).ZT(), 1, "unitimpulse(n)")
        self.assertEqual(3 * unitimpulse(n).ZT(), 3, "unitimpulse(n)")        
        self.assertEqual(unitimpulse(n - 1).ZT(), 1 / z, "unitimpulse(n - 1)")
        self.assertEqual(unitimpulse(n - 2).ZT(), 1 / z**2, "unitimpulse(n - 2)")
        self.assertEqual(3 * unitimpulse(n - 2).ZT(), 3 / z**2, "3 * unitimpulse(n - 2)")
        self.assertEqual(3 * unitimpulse(n + 2).ZT(), 3 * z**2, "3 * unitimpulse(n + 2)")                
        self.assertEqual(expr('v(n)').ZT(), expr('V(z)'), "v(n)")
        self.assertEqual(expr('v(n / 3)').ZT(), expr('V(z**3)'), "v(n/3)")
        self.assertEqual(expr('3 * v(n / 3)').ZT(), expr('3 * V(z**3)'), "3 * v(n/3)")
        self.assertEqual(expr('v(n-3)').ZT(), expr('V(z) / z**3'), "v(n - 3)")
        self.assertEqual(nexpr('x(n)').ZT(), zexpr('X(z)'), "x(n)")

        self.assertEqual((0.1**n).ZT(), 10 * z / (10 * z - 1), "0.1**n")
        self.assertEqual(((0.1**n) * u(n)).ZT(), 10 * z / (10 * z - 1), "0.1**n")        
                         

    def test_inverse_ztransform(self):

        self.assertEqual(zexpr(1).IZT(causal=True), unitimpulse(n), "1")
        self.assertEqual((z**-1).IZT(causal=True), unitimpulse(n-1), "z**-1")
        self.assertEqual((z**-2).IZT(causal=True), unitimpulse(n-2), "z**-2")
        self.assertEqual(zexpr('1 / (1 - a * z ** -1)').IZT(causal=True), nexpr('a**n * u(n)'), "1 / (1 - a * z)")                        
        self.assertEqual(zexpr('X(z)').IZT(causal=True), nexpr('x(n)'), "X(z)")


    def test_misc(self):

        x = 3 * unitimpulse(n - 2)
        y1 = x.differentiate()
        y2 = x.ZT().ndifferentiate().IZT()
        self.assertEqual(y1, y2, "differentiate")

        x = 3 * unitimpulse(n - 2)
        y1 = x.integrate()
        y2 = x.ZT().nintegrate().IZT()
        # Need to teach that these are equal
        #self.assertEqual(y1, y2, "integrate")                
        
    def test_convolution(self):

        a = expr('a(n)')
        b = expr('b(n)')        
        c = (a.ZT() * b.ZT()).IZT(causal=True)
        d = expr('Sum(a(-m + n)*b(m), (m, 0, n))')
        
        self.assertEqual(c, d, "convolution")
