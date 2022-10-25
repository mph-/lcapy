from lcapy import *
from lcapy.jomegaexpr import AngularFrequencyResponseDomainExpression
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

    def test_laplace_jomega_transform(self):

        H = 1 / s
        P = H(jw)

        self.assertEqual(type(P), AngularFrequencyResponseDomainExpression, 'H(jw)')
