from lcapy import *
from lcapy.discretetime import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_DFT(self):

        self.assertEqual(nexpr('delta(n)').DFT(), 1, "delta(n)")
        self.assertEqual(nexpr('2 * delta(n)').DFT(), 2, "2 * delta(n)")

        self.assertEqual(nexpr('1').DFT(), kexpr('N * delta(k)'), "1")
        self.assertEqual(nexpr('2').DFT(), kexpr('2 * N * delta(k)'), "2")

        self.assertEqual(nexpr('x(n)').DFT(), kexpr('X(k)'), "x(n)")
        self.assertEqual(nexpr('2 * x(n)').DFT(), kexpr('2 * X(k)'), "2 * x(n)")        

        self.assertEqual(nexpr('delta(n - 1)').DFT(),
                         kexpr('exp(-j * 2 * pi * k / N)'), "delta(n - 1)")

        self.assertEqual(nexpr('2 * delta(n - 1)').DFT(),
                         kexpr('2 * exp(-j * 2 * pi * k / N)'), "2 * delta(n - 1)")

        self.assertEqual(nexpr('delta(n - 2)').DFT(),
                         kexpr('exp(-j * 2 * pi * 2 * k / N)'), "delta(n - 2)")

        self.assertEqual(nexpr('2 * exp(-j * 2 * pi * n / N)').DFT(), 
                         kexpr('2 * N * delta(k - 1)'), "2 * exp(-j * 2 * pi * n / N)")


    def test_IDFT(self):

        self.assertEqual(kexpr(1).IDFT(), delta(n), "1")
        self.assertEqual(delta(k).IDFT(), nexpr('1 / N'), "delta(k)")                
