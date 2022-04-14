from lcapy import *
import sympy as sym
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_DTFT(self):

        x = delta(n)
        X = x.DTFT(images=0)
        self.assertEqual(X, 1, "delta(n)")
        self.assertEqual(X.IDTFT(), x, "delta(n)")

        y = 2 * delta(n)
        Y = y.DTFT(images=0)
        self.assertEqual(Y, 2, "2 * delta(n)")
        self.assertEqual(Y.IDTFT(), y, "2 * delta(n)")

        Xoo = x.DTFT(images=oo)
        Yoo = y.DTFT(images=oo)
        Zoo = Xoo + Yoo

        self.assertEqual(Zoo.remove_images(), X + Y, "remove_images")
        self.assertEqual(Zoo.IDTFT(), x + y, "sum IDTFT")

        x = 0 * n + 2
        self.assertEqual(x.DTFT(images=0), 2 * delta(f) / dt, "2")
        self.assertEqual(x.DTFT(images=0).IDTFT(), x, "2")

        x = nexpr('x(n)')
        self.assertEqual(x.DTFT(images=0), fexpr('X(f)'), "x(n)")
        self.assertEqual(x.DTFT(images=0).IDTFT(), x, "x(n)")

        x = nexpr('2 * x(n)')
        self.assertEqual(x.DTFT(images=0), fexpr('2 * X(f)'), "2 * x(n)")
        self.assertEqual(x.DTFT(images=0).IDTFT(), x, "2 * x(n)")

        x = nexpr('x(2 * n)')
        self.assertEqual(x.DTFT(images=0), fexpr('X(f / 2) / 2'), "x(2 * n)")
        self.assertEqual(x.DTFT(images=0).IDTFT(), x, "x(2 * n)")

        x = sign(n)
        X = x.DTFT(Omega)
        X1 = 2 * exp(j * Omega) / (exp(j * Omega) - 1)

        self.assertEqual(X, X1, "sign(n).DTFT(Omega)")

        m = symbol('m', integer=True)
        x = nexpr('x(n)')
        X = x.DTFT(f)
        # Produces infinite recursion on replacement
        #Y = Sum(fexpr('X(f)'), (m, -oo, oo)).replace(f, f - m / dt)

        m = symbol('m', integer=True)
        # Produces infinite recursion on comparison
        #self.assertEqual(X, Y, "x(n).DTFT(F)")
