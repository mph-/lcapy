from lcapy import *
from lcapy.expr import TimeDomainExpression
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_hilbert(self):

        self.assertEqual(DiracDelta(t).HT(), 1 / (pi * t), "delta(t)")
        self.assertEqual((1 / t).HT(), pi * DiracDelta(t), "1 / t")
        self.assertEqual(cos(t).HT(), sin(t), "cos(t)")
        self.assertEqual(sin(t).HT(), -cos(t), "cos(t)")
        self.assertEqual(exp(j * t).HT(), -j * exp(j * t), "exp(j * t)")
        self.assertEqual(exp(-j * t).HT(), j * exp(-j * t), "exp(-j * t)")

    def test_inverse_hilbert(self):

        self.assertEqual((1 / t).IHT(), pi * DiracDelta(t), "-1 / t")
