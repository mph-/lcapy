from lcapy import *
import unittest

class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy"""

    def assertEqual2(self, ans1, ans2, comment):

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            pprint(ans1)
            pprint(ans2)
            raise AssertionError(e)

    def test_approximate_exp(self):

        H = exp(-s)
        Ha = H.approximate_exp(order=1)
        Hb = H.approximate_exp(order=2, numer_order=1)

        self.assertEqual2(Ha, (2 - s)/ (2 + s), 'approximate exp Pade order 1')
        self.assertEqual2(Hb, (6 - 2 * s)/ (6 + 4 *s + s**2), 'approximate exp Pade order 1,2')
