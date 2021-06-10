from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_DTFT(self):

        self.assertEqual(nexpr('delta(n)').DTFT(), 1, "delta(n)")
        self.assertEqual(nexpr('2 * delta(n)').DTFT(), 2, "2 * delta(n)")

        self.assertEqual(nexpr('2').DTFT(), fexpr('2 * delta(f)'), "2")

        self.assertEqual(nexpr('x(n)').DTFT(), fexpr('X(f)'), "x(n)")
        self.assertEqual(nexpr('2 * x(n)').DTFT(), fexpr('2 * X(f)'), "2 * x(n)")
        self.assertEqual(nexpr('x(2 * n)').DTFT(), fexpr('X(f / 2) / 2'), "x(2 * n)")        
