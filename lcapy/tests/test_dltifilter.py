from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_DLTIFilter(self):

        fil = DLTIFilter([1] * 2, [1])

        self.assertEqual(fil.impulse_response(), UnitImpulse(
            n) + UnitImpulse(n - 1), "impulse_response")
        self.assertEqual(fil.transfer_function(),
                         (z + 1) / z, "transfer_function")
        self.assertEqual(fil.difference_equation(), DifferenceEquation(
            expr('y(n)'), expr('x(n) + x(n - 1)')), "difference_equation")
