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
        
