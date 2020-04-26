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
        
        self.assertEqual(b, 6 * pi * f, 'transform omega -> f')
        self.assertEqual(c, 6 * pi * f, 'explicit transform omega -> f')

    def test_f_to_omega(self):

        a = 2 * pi * f
        b = a(omega)
        c = a.transform(omega)
        
        self.assertEqual(b, omega, 'transform f -> omega')
        self.assertEqual(c, omega, 'explicit transform f -> omega')

    def test_f_to_f(self):

        a = 2 * pi * f
        b = a.transform(f)
        
        self.assertEqual(a, b, 'explicit transform f -> f')

    def test_omega_to_omega(self):

        a = 2 * omega
        b = a.transform(omega)
        
        self.assertEqual(a, b, 'explicit transform omega -> omega')        

    def test_noisef_to_noiseomega(self):

        a = Vnoisy(omega)
        b = a(omega)
        c = a(f)
        d = c(omega)

        self.assertEqual(a, b, 'noisyomega -> noisyomega')
        self.assertEqual(a, d, 'noisyomega -> noisyf -> noisyomega')
        
    def test_s_to_jomega(self):

        a = expr(2 * s, causal=True)
        b = a.transform(jomega)

        self.assertEqual(b, 2 * jomega, 'explicit transform s -> jomega')

        a = 2 * s
        b = a.transform(jomega, causal=True)        
        self.assertEqual(b, 2 * jomega, 'explicit transform s -> jomega')

    def test_s_to_omega(self):

        a = expr(2 * s, causal=True)
        b = a.transform(omega)

        self.assertEqual(b, 2 * j * omega, 'explicit transform s -> omega')

        a = 2 * s
        b = a.transform(omega, causal=True)        
        self.assertEqual(b, 2 * j * omega, 'explicit transform s -> omega')

    def test_s_to_f(self):

        a = expr(2 * s, causal=True)
        b = a.transform(f)

        self.assertEqual(b, j * 4 * pi * f, 'explicit transform s -> f')

        a = 2 * s
        b = a.transform(f, causal=True)        
        self.assertEqual(b, j * 4 * pi * f, 'explicit transform s -> f')
        
        
