from lcapy import *
from lcapy.laplace import inverse_laplace_ratfun
from lcapy.sexpr import LaplaceDomainVoltage
from lcapy.texpr import TimeDomainExpression, TimeDomainVoltage
import unittest
import sympy as sym


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_laplace(self):

        a = 0 * t + 1
        r1, r2 = inverse_laplace_ratfun(a.expr, s.var, t.var)
        self.assertEqual(r1, delta(t).expr, "inverse_laplace_ratfun")
        self.assertEqual(r2, 0, "inverse_laplace_ratfun")        

        a = 1 / (s + 2)
        r1, r2 = inverse_laplace_ratfun(a.expr, s.var, t.var)
        self.assertEqual(r1, 0, "inverse_laplace_ratfun")        
        self.assertEqual(r2, exp(-2 * t).expr, "inverse_laplace_ratfun")

        self.assertEqual(Heaviside(t).laplace(), 1 / s, "Heaviside(t)")
        self.assertEqual(DiracDelta(t).laplace(), 1, "DiracDelta(t)")
        self.assertEqual(TimeDomainVoltage('x(t)').laplace(), LaplaceDomainVoltage('X(s)'), "x(t)")
        self.assertEqual(TimeDomainVoltage('5 * x(t)').laplace(), LaplaceDomainVoltage('5 * X(s)'), "5 * x(t)")

        v = expr('R0 * exp(-alpha * t) * i(t)') 
        V = expr('R0 * I(s + alpha)')

        self.assertEqual(v.laplace(), V, "R0 * exp(-alpha * t) * i(t)")
        
        
    def test_inverse_laplace(self):

        self.assertEqual((1 / s).inverse_laplace(causal=True), Heaviside(t),
                         "1 / s")
        self.assertEqual((s * 0 + 1).inverse_laplace(causal=True), DiracDelta(t),
                         "1")
        self.assertEqual((s * 0 + 10).inverse_laplace(causal=True), 10
                         * DiracDelta(t), "0")
        self.assertEqual(LaplaceDomainVoltage('V(s)').inverse_laplace(causal=True),
                         TimeDomainVoltage('v(t)'), "V(s)")
        self.assertEqual(LaplaceDomainVoltage('10 * V(s)').inverse_laplace(causal=True),
                         TimeDomainVoltage('10 * v(t)'), "V(s)")
        self.assertEqual(LaplaceDomainVoltage('10 * V(s) * exp(-5 * s)').inverse_laplace(causal=True), TimeDomainVoltage('10 * v(t - 5)'), "10 * V(s) * exp(-5 * s)")
        self.assertEqual(TimeDomainVoltage('v(t)').laplace().inverse_laplace(causal=True),
                         TimeDomainVoltage('v(t)'), "v(t)")
        self.assertEqual(expr('1/(s+a)').inverse_laplace(causal=True), expr('exp(-a * t) * u(t)'), "1/(s+a)")
        self.assertEqual(expr('1/(s**2)').inverse_laplace(causal=True), expr('t * u(t)'), "1/(s**2)")        
        self.assertEqual(expr('s/(s+a)').inverse_laplace(causal=True), expr('-a * exp(-a * t) * u(t) + delta(t)'), "s/(s+a)")
        self.assertEqual(expr('s/(s**2+a**2)').inverse_laplace(causal=True), expr('cos(a * t) * u(t)'), "s/(s**2+a**2)")
        self.assertEqual(expr('a/(s**2+a**2)').inverse_laplace(causal=True), expr('sin(a * t) * u(t)'), "a/(s**2+a**2)")                                                                           

    def test_damped_sin(self):

        H1 = 2 / (2 * s ** 2 + 5 * s + 6)
        H2 = H1 * s
        H3 = H1 * s * s

        self.assertEqual(H1(t, damped_sin=True)(s), H1, "damped sin1")
        self.assertEqual(H2(t, damped_sin=True)(s), H2, "damped sin2")
        self.assertEqual(H3(t, damped_sin=True)(s), H3, "damped sin3")
        self.assertEqual((H1 + H2)(t, damped_sin=True)(s), H1 + H2, "damped sin1, 2")
        self.assertEqual((H1 + H3)(t, damped_sin=True)(s), H1 + H3, "damped sin1, 3")
        self.assertEqual((H1 + H2 + H3)(t, damped_sin=True)(s), H1 + H2 + H3, "damped sin1, 2, 3")                


    def test_derivative_undef(self):

        H = s * 'I(s)'
        h = H(t)
        H2 = h(s)

        self.assertEqual(H, H2, "derivative of undef")

        H = s**2 * 'I(s)'
        h = H(t)
        H2 = h(s)

        self.assertEqual(H, H2, "second derivative of undef")                

    def test_laplace_convolution(self):

        a = expr('Integral(3 * x(t - tau) * y(tau), (tau, -oo, oo))')
        self.assertEqual(a(s), expr('3 * X(s) * Y(s)'), "3 * X * Y")
        a = expr('Integral(3 * x(tau) * y(t - tau), (tau, -oo, oo))')
        self.assertEqual(a(s), expr('3 * X(s) * Y(s)'), "3 * X * Y")
        a = expr('Integral(3 * x(t - tau) * Heaviside(tau), (tau, -oo, oo))')
        self.assertEqual(a(s), expr('3 * X(s) / s'), "3 * X / s")
        a = expr('Integral(3 * x(tau) * Heaviside(t - tau), (tau, -oo, oo))')
        self.assertEqual(a(s), expr('3 * X(s) / s'), "3 * X / s")

        a = expr('Integral(3 * x(t - tau), (tau, 0, oo))')
        self.assertEqual(a(s), expr('3 * X(s) / s'), "3 * X / s")
        
 
