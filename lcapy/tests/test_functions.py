from lcapy import *
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

        self.assertEqual(rect(n).rewrite(), us(n + 1 / 2) - us(n - 1 / 2), "rect(n)")
        

        
