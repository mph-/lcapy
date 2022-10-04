from lcapy import *
from lcapy.sexpr import LaplaceDomainVoltage
from lcapy.texpr import TimeDomainExpression, TimeDomainVoltage
import unittest
import sympy as sym


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_laplace(self):

        self.assertEqual(Heaviside(t).laplace(), 1 / s, "Heaviside(t)")
        self.assertEqual(DiracDelta(t).laplace(), 1, "DiracDelta(t)")
        self.assertEqual(TimeDomainVoltage('x(t)').laplace(),
                         LaplaceDomainVoltage('X(s)'), "x(t)")
        self.assertEqual(TimeDomainVoltage('5 * x(t)').laplace(),
                         LaplaceDomainVoltage('5 * X(s)'), "5 * x(t)")

        v = expr('R0 * exp(-alpha * t) * i(t)')
        V = expr('R0 * I(s + alpha)')

        self.assertEqual(v.laplace(), V, "R0 * exp(-alpha * t) * i(t)")

        self.assertEqual(rect(t - 0.5).laplace(), (1 - exp(-s)) / s, "rect(t)")
        self.assertEqual(rect(t / 2 - 0.5).laplace(),
                         (1 - exp(-2 * s)) / s, "rect(t / 2)")
        self.assertEqual(ramp(t).laplace(), 1 / s**2, "ramp(t)")
        self.assertEqual(rampstep(t).laplace(),
                         (1 - exp(-s)) / s**2, "rampstep(t)")

    def test_inverse_laplace(self):

        self.assertEqual((1 / s).inverse_laplace(causal=True), Heaviside(t),
                         "1 / s")
        self.assertEqual((s * 0 + 1).inverse_laplace(causal=True), DiracDelta(t),
                         "1")
        self.assertEqual((s).inverse_laplace(causal=True), DiracDelta(t, 1),
                         "s")
        self.assertEqual((s * 0 + 10).inverse_laplace(causal=True), 10
                         * DiracDelta(t), "10")
        self.assertEqual(LaplaceDomainVoltage('V(s)').inverse_laplace(causal=True),
                         TimeDomainVoltage('v(t)'), "V(s)")
        self.assertEqual(LaplaceDomainVoltage('10 * V(s)').inverse_laplace(causal=True),
                         TimeDomainVoltage('10 * v(t)'), "V(s)")
        self.assertEqual(LaplaceDomainVoltage('10 * V(s) * exp(-5 * s)').inverse_laplace(
            causal=True), TimeDomainVoltage('10 * v(t - 5)'), "10 * V(s) * exp(-5 * s)")
        self.assertEqual(TimeDomainVoltage('v(t)').laplace().inverse_laplace(causal=True),
                         TimeDomainVoltage('v(t)'), "v(t)")
        self.assertEqual(expr('1/(s+a)').inverse_laplace(causal=True),
                         expr('exp(-a * t) * u(t)'), "1/(s+a)")
        self.assertEqual(
            expr('1/(s**2)').inverse_laplace(causal=True), expr('t * u(t)'), "1/(s**2)")
        self.assertEqual(expr('1/((s+3)**2)').inverse_laplace(causal=True),
                         expr('t * u(t) * exp(-3 * t)'), "1/((s+3)**2)")
        self.assertEqual(expr('1/(s**3)').inverse_laplace(causal=True),
                         expr('t**2 * u(t) / 2'), "1/(s**3)")
        self.assertEqual(expr('s/(s+a)').inverse_laplace(causal=True),
                         expr('-a * exp(-a * t) * u(t) + delta(t)'), "s/(s+a)")
        self.assertEqual(expr('s/(s**2+a**2)').inverse_laplace(causal=True),
                         expr('cos(a * t) * u(t)'), "s/(s**2+a**2)")
        self.assertEqual(expr('a/(s**2+a**2)').inverse_laplace(causal=True),
                         expr('sin(a * t) * u(t)'), "a/(s**2+a**2)")
        self.assertEqual(expr('(1/s - exp(-s)/s)/s').ILT(causal=True),
                         expr('(t*u(t) + (1 - t)*u(t - 1))'), "(1/s - exp(-s)/s)/s")
        self.assertEqual(expr('(1/s - exp(-s)/s)/s').expand().ILT(causal=True),
                         expr('(t*u(t) + (1 - t)*u(t - 1))'), "(1/s - exp(-s)/s)/s")

        self.assertEqual(expr('s*(1/s - exp(-5*s)/s)/(s + 1)').ILT(causal=True),
                         expr('-exp(5 - t)*u(t - 5) + exp(-t)*u(t)'),
                         "'s*(1/s - exp(-5*s)/s)/(s + 1)")

        self.assertEqual(expr('s*(1/s - exp(-5*s)/s)/(s + 1)').expand().ILT(causal=True),
                         expr('-exp(5 - t)*u(t - 5) + exp(-t)*u(t)'),
                         "'s*(1/s - exp(-5*s)/s)/(s + 1)")

        A4 = (1/s - exp(-5*s)/s)/(s**2 + 707*s/500 + 1)
        self.assertEqual(A4(t)(s), A4, str(A4))
        self.assertEqual(expr('a/(s+a)', causal=True)(t),
                         expr('a*exp(-a * t)*u(t)'), "s/(s+a)")

        A5 = (s + 1) * (s + 2) * expr('X(s)')
        self.assertEqual(A5(t)(s), A5, str(A5))

    def test_damped_sin(self):

        H1 = 2 / (2 * s ** 2 + 5 * s + 6)
        H2 = H1 * s
        H3 = H1 * s * s

        self.assertEqual(H1(t, damped_sin=True)(s), H1, "damped sin1")
        self.assertEqual(H2(t, damped_sin=True)(s), H2, "damped sin2")
        self.assertEqual(H3(t, damped_sin=True)(s), H3, "damped sin3")
        self.assertEqual((H1 + H2)(t, damped_sin=True)
                         (s), H1 + H2, "damped sin1, 2")
        self.assertEqual((H1 + H3)(t, damped_sin=True)
                         (s), H1 + H3, "damped sin1, 3")
        self.assertEqual((H1 + H2 + H3)(t, damped_sin=True)
                         (s), H1 + H2 + H3, "damped sin1, 2, 3")

    def test_derivative_undef(self):

        H = s * 'I(s)' - 'i(0)'
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

        a = expr('Integral(3 * x(tau), (tau, 0, t))')
        self.assertEqual(a(s), expr('3 * X(s) / s'), "3 * X / s")

        a = expr('x(t)').convolve(expr('g(t)'))
        self.assertEqual(a(s), expr('X(s) * G(s)'), "X(s) * G(s)")

        a = expr('x(t)').convolve(expr('g(t)'), commutate=True)
        self.assertEqual(a(s), expr('X(s) * G(s)'), "X(s) * G(s)")

        a = expr('x(t)').convolve(expr('g(t)', causal=True))
        self.assertEqual(a(s), expr('X(s) * G(s)'), "X(s) * G(s)")

        a = expr('x(t)').convolve(expr('g(t)', causal=True), commutate=True)
        self.assertEqual(a(s), expr('X(s) * G(s)'), "X(s) * G(s)")

        a = expr('x(t)', causal=True).convolve(expr('g(t)'))
        self.assertEqual(a(s), expr('X(s) * G(s)'), "X(s) * G(s)")

        a = expr('x(t)', causal=True).convolve(expr('g(t)', commutate=True))
        self.assertEqual(a(s), expr('X(s) * G(s)'), "X(s) * G(s)")

        a = expr('x(t)', causal=True).convolve(expr('g(t)', causal=True))
        self.assertEqual(a(s), expr('X(s) * G(s)'), "X(s) * G(s)")

        a = expr('x(t)', causal=True).convolve(expr('g(t)', causal=True,
                                                    commutate=True))
        self.assertEqual(a(s), expr('X(s) * G(s)'), "X(s) * G(s)")
