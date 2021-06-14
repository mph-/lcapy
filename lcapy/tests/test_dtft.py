from lcapy import *
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
        self.assertEqual(X, 2 / (1 - exp(-j * Omega)), "sign(n).DTFT(Omega)")
