from lcapy import *
from lcapy.laplace import inverse_laplace_ratfun
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_laplace(self):

        a = 0 * t + 1
        r1, r2 = inverse_laplace_ratfun(a.expr, s.var, t.var)
        self.assertEqual(r1, delta(t), "inverse_laplace_ratfun")
        self.assertEqual(r2, 0, "inverse_laplace_ratfun")        

        a = 1 / (s + 2)
        r1, r2 = inverse_laplace_ratfun(a.expr, s.var, t.var)
        self.assertEqual(r1, 0, "inverse_laplace_ratfun")        
        self.assertEqual(r2, exp(-2 * t).expr, "inverse_laplace_ratfun")

        self.assertEqual(Heaviside(t).laplace(), 1 / s, "Heaviside(t)")
        self.assertEqual(DiracDelta(t).laplace(), 1, "DiracDelta(t)")
        self.assertEqual(Vt('x(t)').laplace(), Vs('X(s)'), "x(t)")
        self.assertEqual(Vt('5 * x(t)').laplace(), Vs('5 * X(s)'), "5 * x(t)")

    def test_inverse_laplace(self):

        self.assertEqual((1 / s).inverse_laplace(causal=True), Heaviside(t),
                         "1 / s")
        self.assertEqual((s * 0 + 1).inverse_laplace(causal=True), DiracDelta(t),
                         "1")
        self.assertEqual((s * 0 + 10).inverse_laplace(causal=True), 10
                         * DiracDelta(t), "0")
        self.assertEqual(Vs('V(s)').inverse_laplace(causal=True),
                         Vt('v(t)'), "V(s)")
        self.assertEqual(Vs('10 * V(s)').inverse_laplace(causal=True),
                         Vt('10 * v(t)'), "V(s)")
        self.assertEqual(Vs('10 * V(s) * exp(-5 * s)').inverse_laplace(causal=True), Vt('10 * v(t - 5)'), "10 * V(s) * exp(-5 * s)")
        self.assertEqual(Vt('v(t)').laplace().inverse_laplace(causal=True),
                         Vt('v(t)'), "v(t)")
                         
