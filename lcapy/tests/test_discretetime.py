from lcapy import *
from lcapy.discretetime import *
import unittest
import sympy as sym


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_transforms(self):

        a = ui(n)
        self.assertEqual(a(z), 1, "ui(n)")
        self.assertEqual(a(k), 1, "ui(k)")
        self.assertEqual(a(k)(z), 1, "ui(k)(z)")
        self.assertEqual(a(k)(n), a, "ui(k)(n)")
        self.assertEqual(a(z)(n), a, "ui(z)(n)")
        self.assertEqual(a(z)(k), 1, "ui(z)(k)")

    def test_as_AB(self):

        H = (z - 2) / (z - 3)
        A, B = H.as_AB()

        self.assertEqual(A, 1 - 3 / z, "A")
        self.assertEqual(B, 1 - 2 / z, "A")

    def test_sequence(self):

        x = seq((1, 2, 3))
        e = delta(n) + 2 * delta(n - 1) + 3 * delta(n - 2)

        y = seq((1, 2, 3, 0, 0))

        h = seq((1, 2))

        self.assertEqual(x.as_impulses(), e, "as_impulses")
        self.assertEqual(x.zeropad(2), y, "zeropad")
        self.assertEqual(x.prune(), x, "prune")
        self.assertEqual(y.extent, 3, "extent")
        self.assertEqual(y.origin, 0, "origin")
        self.assertEqual(x.convolve(h), seq((1, 4, 7, 6)), "convolve")
        self.assertEqual(h.expr.ZT(), (z + 2) / z, "ZT")
        self.assertEqual(
            h.latex(), r'\left\{\underline{1}, 2\right\}', "latex")
        self.assertEqual(seq('1, 2, 3'), x, "seq")
        self.assertEqual(e.seq(), x, "seq")

        self.assertEqual(h.n, [0, 1], "seq.n")
        self.assertEqual(h.expr(z), (z + 2) / z, "h(z)")

        self.assertEqual(h.delay(1), seq((1, 2), origin=-1), "delay")
        self.assertEqual(h.delay(1), h >> 1, ">>")
        self.assertEqual(h.delay(-1), h << 1, "<<")

        q = seq('{1, _2, 3, 4}')
        self.assertEqual(q.extent, 4, "extent")
        self.assertEqual(q.origin, 1, "origin")
        self.assertEqual(q.n, [-1, 0, 1, 2], "seq.n")
        self.assertEqual(q[0], 2, "seq[0]")
        self.assertEqual(
            q.latex(), r'\left\{1, \underline{2}, 3, 4\right\}', "latex")
        q.origin = 2
        self.assertEqual(q.n, [-2, -1, 0, 1], "seq.n with origin 2")

        a = seq((1, 2, 3), origin=3)
        self.assertEqual(a.zeroextend().vals, [1, 2, 3, 0], "zeroextend")

        x = seq((1, 2, 3, 4))
        X = x.DFT()
        x2 = X.IDFT()
        self.assertEqual(x, x2, "DFT/IDFT")

        self.assertEqual(h.ZT(), zseq((1, 2 / z)), "ZT")
        self.assertEqual(h.ZT().IZT(), h, "ZT/IZT")

    def test_zexpr(self):

        H = z / (z - 1)

        de = H.difference_equation().separate()
        rhs = expr('x(n)')
        lhs = expr('y(n) - y(n-1)')

        self.assertEqual(de.lhs, lhs, "DE lhs")
        self.assertEqual(de.rhs, rhs, "DE rhs")

    def test_nexpr(self):

        x = seq('{1, _2, 3, 4}').as_impulses()

        y1 = seq('{0, 1, 2, 3, 4}').as_impulses()
        y2 = x.delay(2)

        self.assertEqual(y1, y2, "delay(2)")
        self.assertEqual(x.first_index(), -1, "first_index")
        self.assertEqual(x.last_index(), 2, "last_index")

        x1 = exp('n')
        x2 = exp(n)
        self.assertEqual(x1, x2, "expr")

    def test_undef_solve(self):
        e = expr('Eq(y(n), a * y(n - 1) + (1 - a) * x(n))')
        E = e(z)
        Y = E.solve('Y')[0]

        Yans = expr('z * (a - 1) * X(z) / (a - z)')

        self.assertEqual(Y, Yans, "solve with undefs")

    def test_dtfilter(self):

        a = symbol('a')
        lpf = DLTIFilter((1 - a, ), (1, -a))

        y = lpf.response(1, [0], (0, 2))

        self.assertEqual(y[0], 1 - a, "IIR impulse response @0")
        self.assertEqual(y[1], a * (1 - a), "IIR impulse response @1")
        self.assertEqual(y[2], a**2 * (1 - a), "IIR impulse response @2")

        y = lpf.response(0, [1], (0, 2))

        self.assertEqual(y[0], a, "transient response @0")
        self.assertEqual(y[1], a**2, "transient response @1")

    def test_dt_assumptions(self):

        self.assertEqual((1 + 1 / z).is_dc, False, "is_dc")
