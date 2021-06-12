from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_DTFT(self):

        x = delta(n)
        self.assertEqual(x.DTFT(), 1, "delta(n)")
        self.assertEqual(x.DTFT().IDTFT(), x, "delta(n)")        

        x = 2 * delta(n)
        self.assertEqual(x.DTFT(), 2, "2 * delta(n)")
        self.assertEqual(x.DTFT().IDTFT(), x, "2 * delta(n)")        

        x = 0 * n + 2
        self.assertEqual(x.DTFT(), 2 * delta(f), "2")
        self.assertEqual(x.DTFT().IDTFT(), x, "2")        

        x = nexpr('x(n)')
        self.assertEqual(x.DTFT(), fexpr('X(f)'), "x(n)")
        self.assertEqual(x.DTFT().IDTFT(), x, "x(n)")

        x = nexpr('2 * x(n)')        
        self.assertEqual(x.DTFT(), fexpr('2 * X(f)'), "2 * x(n)")
        self.assertEqual(x.DTFT().IDTFT(), x, "2 * x(n)")        

        x = nexpr('x(2 * n)')                
        self.assertEqual(x.DTFT(), fexpr('X(f / 2) / 2'), "x(2 * n)")
        self.assertEqual(x.DTFT().IDTFT(), x, "x(2 * n)")        

        
