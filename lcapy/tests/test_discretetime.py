from lcapy import *
from lcapy.discretetime import *
import unittest
import sympy as sym


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_transform(self):

        a = ui(n)

        self.assertEqual(a(z), 1, "ui(n)")

    def test_decompose_AB(self):

        H = (z -2) / (z - 3)
        A, B = H.decompose_AB()

        self.assertEqual(A, 1 - 3 / z, "A")
        self.assertEqual(B, 1 - 2 / z, "A")                

    def test_sequence(self):

        x = seq((1, 2, 3))
        e = delta(n) + 2 * delta(n - 1) + + 3 * delta(n - 2)

        y = seq((1, 2, 3, 0, 0))

        h = seq((1, 2))                
        
        self.assertEqual(x.as_impulses(), e, "as_impulses")
        self.assertEqual(x.zeropad(2), y, "zeropad")
        self.assertEqual(x.prune(), x, "prune")
        self.assertEqual(y.extent(), 3, "extent")
        self.assertEqual(x.convolve(h), seq((1, 4, 7, 6)), "convolve")
        self.assertEqual(h.ZT(), (z + 2) / z, "ZT")
        self.assertEqual(h.latex(), r'\left\{\underline{1}, 2\right\}', "latex")                
        self.assertEqual(seq('1, 2, 3'), x, "seq")
        self.assertEqual(seq('{1, _2, 3}').latex(), r'\left\{1, \underline{2}, 3\right\}', "latex")
        
        
