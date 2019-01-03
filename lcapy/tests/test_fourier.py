from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_fourier(self):

        self.assertEqual(Heaviside(t).fourier(), 1 / s, "Heaviside(t)")
        self.assertEqual(DiracDelta(t).fourier(), 1, "DiracDelta(t)")
        self.assertEqual(Vt('x(t)').fourier(), Vf('X(f)'), "x(t)")
        self.assertEqual(Vt('5 * x(t)').fourier(), Vf('5 * X(f)'), "5 * x(t)")

    def test_inverse_fourier(self):

        self.assertEqual((f * 0 + 1).inverse_fourier(), DiracDelta(t),
                         "1")
        self.assertEqual((f * 0 + 10).inverse_fourier(), 10
                         * DiracDelta(t), "0")
        self.assertEqual(Vf('V(f)').inverse_fourier(), Vt('v(t)'), "V(f)")
        self.assertEqual(Vf('10 * V(f)').inverse_fourier(),
                         Vt('10 * v(t)'), "V(f)")
        self.assertEqual(Vt('v(t)').fourier().inverse_fourier(),
                         Vt('v(t)'), "v(t)")
                         
