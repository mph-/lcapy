from lcapy import R, C, L, V, I, v, exp, Heaviside, Vac, Vdc
from lcapy.core import Zs, s, t
import unittest
import sympy as sym


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def assertEqual2(self, ans1, ans2, comment):

        ans1 = ans1.canonical()
        ans2 = ans2.canonical()

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            ans1.pprint()
            ans2.pprint()
            raise AssertionError(e)

    def test_VR_ac1(self):
        """Lcapy: check VR ac network

        """

        a = Vac(4) | R(2)

        self.assertEqual(a.is_dc, False, "DC incorrect")
        self.assertEqual(a.is_ac, True, "AC incorrect")
        self.assertEqual(a.i, 0, "Open circuit terminal current not zero")


    def test_VR_dc1(self):
        """Lcapy: check VR dc network

        """

        a = Vdc(4) | R(2)

        self.assertEqual(a.is_ac, False, "AC incorrect")
        self.assertEqual(a.is_dc, True, "DC incorrect")
        self.assertEqual(a.i, 0, "Open circuit terminal current not zero")

    def test_VC_dc1(self):
        """Lcapy: check VC series dc network

        """

        a = Vdc(4) + C(2)

        self.assertEqual(a.is_ac, False, "AC incorrect")
        self.assertEqual(a.is_dc, True, "DC incorrect")
        self.assertEqual(a.i, 0, "Open circuit terminal current not zero")
        self.assertEqual(a.voc, 4, "Open circuit terminal voltage incorrect")
