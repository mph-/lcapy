from lcapy import *
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


    def test_VR_dc1(self):
        """Lcapy: check VR dc network

        """

        a = Vdc(4) | R(2)

        self.assertEqual(a.initial_value_problem, False, "initial_value_problem incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")
        self.assertEqual(a.is_dc, True, "DC incorrect")


    def test_VC_dc1(self):
        """Lcapy: check VC dc network

        """

        a = Vdc(4) + C(2)

        self.assertEqual(a.initial_value_problem, False, "initial_value_problem incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")
        self.assertEqual(a.is_dc, True, "DC incorrect")


    def test_VC_dc2(self):
        """Lcapy: check VC dc network

        """

        a = Vdc(4) + C(2, 0)

        self.assertEqual(a.initial_value_problem, True, "initial_value_problem incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")
        self.assertEqual(a.is_dc, False, "DC incorrect")

    def test_thevenin_ac(self):
        """Lcapy: check ac Thevenin conversion

        """
        a = (Vac('1') + C(2)) | R(3)
        
        self.assertEqual(a.norton().isc, It(2 * omega * sin(omega * t)),
                         "Isc incorrect")
        self.assertEqual(a.thevenin().isc,  It(2 * omega * sin(omega * t)),
                         "Isc incorrect")
        
    def test_superposition(self):
        """Lcapy: check network superposition"""

        a = Vac(40) + Vnoise(20) + Vstep(10) + R(5)
        self.assertEqual2(a.Voc.s, 10 / s, "Voc.s error")
        self.assertEqual2(a.Voc.n, Vn(20), "Voc.n error")
        self.assertEqual2(a.Voc.w, Vphasor(40), "Voc.w error")
        self.assertEqual2(a.Isc.s, 2 / s, "Isc.s error")
        # FIXME, this intermittently fails.
        self.assertEqual2(a.Isc.n, In(4), "Isc.n error")
        self.assertEqual2(a.Isc.w, Iphasor(8), "Isc.w error")        
        
    def test_ivp(self):
        """Lcapy: check network with initial values"""

        a = Vstep(10) + C('C1', 5)
        self.assertEqual(a.is_ivp, True, "IVP fail")
        self.assertEqual2(a.Voc.s, 15 / s, "Voc fail")        
        
