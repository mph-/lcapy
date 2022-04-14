from lcapy import *
from lcapy.sym import tausym
from lcapy.texpr import TimeDomainExpression
from lcapy.omegaexpr import AngularFourierDomainExpression
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_functions(self):

        self.assertEqual(sinc(0), 1, "sinc(0)")
        self.assertEqual(sinc(1), 0, "sinc(1)")
        self.assertEqual(rect(0), 1, "rect(0)")
        self.assertEqual(rect(1), 0, "rect(1)")
        self.assertEqual(ui(1), 0, "ui(1)")
        self.assertEqual(ui(0), 1, "ui(0)")
        self.assertEqual(us(1), 1, "us(1)")
        self.assertEqual(us(0), 1, "us(0)")
        self.assertEqual(us(-1), 0, "us(-1)")

        self.assertEqual(rect(n).rewrite(), us(
            n + 1 / 2) - us(n - 1 / 2), "rect(n)")

        r = Integral(u(t - tausym), (tausym, 0, t))
        e = u(t).convolve(u(t)).simplify()

        self.assertEqual(e, r, "simplify integral with u(t)")

        f0 = symbol('f_0')
        c = cos(2 * pi * f0 * t)
        self.assertEqual(type(c), TimeDomainExpression, "cos")

        a = atan2(1, omega)
        self.assertEqual(type(a), AngularFourierDomainExpression, "atan2")

        self.assertEqual(sign(-n).evaluate(0), 1, "sign(-n)")
        self.assertEqual(sign(-n).evaluate(1), -1, "sign(-n)")

    def test_heaviside(self):

        self.assertEqual((u(t) * 3 * u(t)).simplify(),
                         3 * u(t), "u(t) * 3 * u(t)")

    def test_dirac_delta_simplify(self):

        self.assertEqual(expr('delta(t) * exp(-3 * t)').simplify(),
                         expr('delta(t)'),
                         "delta(t) * exp(-3 * t)")

        self.assertEqual(expr('exp(-3 * t) * delta(t)').simplify(),
                         expr('delta(t)'),
                         "exp(-3 * t) * delta(t)")

        self.assertEqual((rect(t / 4) * delta(t - 1)).simplify(),
                         delta(t - 1),
                         "rect(t / 4) * delta(t - 1)")

    def test_sinc(self):

        a = sinc(n).evaluate(2)
        self.assertEqual(abs(a) < 1e-16, True, "sinc evaluate incorrect.")

        a = sincn(n).evaluate(2)
        self.assertEqual(abs(a) < 1e-16, True, "sincn evaluate incorrect.")

        a = expr('sinc(n)').evaluate(2)
        self.assertEqual(abs(a) < 1e-16, True, "'sinc' evaluate incorrect.")

        a = expr('sincn(n)').evaluate(2)
        self.assertEqual(abs(a) < 1e-16, True, "'sincn' evaluate incorrect.")
