from lcapy import *
from lcapy.texpr import TimeDomainVoltage, TimeDomainCurrent
from lcapy.noiseomegaexpr import AngularFourierNoiseDomainVoltage, AngularFourierNoiseDomainCurrent
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

        self.assertEqual(a.is_ivp, False, "is_ivp incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")
        self.assertEqual(a.is_dc, True, "DC incorrect")

    def test_VC_dc1(self):
        """Lcapy: check VC dc network

        """

        a = Vdc(4) + C(2)

        self.assertEqual(a.is_ivp, False, "is_ivp incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")
        self.assertEqual(a.is_dc, True, "DC incorrect")

    def test_VC_dc2(self):
        """Lcapy: check VC dc network

        """

        a = Vdc(4) + C(2, 0)

        self.assertEqual(a.is_ivp, True, "is_ivp incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")
        self.assertEqual(a.is_dc, False, "DC incorrect")

    def test_thevenin_ac(self):
        """Lcapy: check ac Thevenin conversion

        """
        a = (Vac('1') + C(2)) | R(3)
        
        self.assertEqual(a.norton().isc, TimeDomainCurrent(-2 * omega0 * sin(omega0 * t)),
                         "Isc incorrect")
        self.assertEqual(a.thevenin().isc, TimeDomainCurrent(-2 * omega0 * sin(omega0 * t)),
                         "Isc incorrect")
        
    def test_superposition(self):
        """Lcapy: check network superposition"""

        a = Vac(40) + Vnoise(20) + Vstep(10) + R(5)
        self.assertEqual(a.Voc.dc, 0, "Voc.dc error")
        self.assertEqual(a.Voc.has_dc, False, "Voc.has_dc error")
        self.assertEqual(a.Voc.is_dc, False, "Voc.is_dc error")
        self.assertEqual(a.Voc.has_ac, True, "Voc.has_ac error")
        self.assertEqual(a.Voc.is_ac, False, "Voc.is_ac error")                                
        self.assertEqual2(a.Voc.s, voltage(10 / s), "Voc.s error")
        self.assertEqual(a.Voc.n.expr, AngularFourierNoiseDomainVoltage(20).expr, "Voc.n error")
        self.assertEqual2(a.Isc.s, current(2 / s), "Isc.s error")
        # FIXME, this intermittently fails.
        self.assertEqual(a.Isc.n.expr, AngularFourierNoiseDomainCurrent(4).expr, "Isc.n error")
        
    def test_ivp(self):
        """Lcapy: check network with initial values"""

        a = Vstep(10) + C('C1', 5)
        self.assertEqual(a.is_ivp, True, "IVP fail")
        self.assertEqual2(a.Voc.s, voltage(15 / s), "Voc fail")

    def test_causal(self):
        """Lcapy: check network is causal"""

        a = Vstep(10) + C('C1')
        self.assertEqual(a.is_causal, True, "causal fail")
        self.assertEqual(a.Isc.is_causal, True, "causal fail")        
        
    def test_YZ(self):

        a = R(1) + V(2) + R(3)
        self.assertEqual(a.Z, impedance(4), "series Z")
        self.assertEqual(a.Y, admittance('1 / 4'), "series Y")

        b = (R(1) + V(2)) | R(3)
        self.assertEqual(b.Z, impedance('3 / 4'), "par Z")
        self.assertEqual(b.Y, admittance('4 / 3'), "par Y")

    def test_attrs(self):

        self.assertEqual(R(1).is_resistor, True, "is_resistor")
        self.assertEqual(G(1).is_conductor, True, "is_conductor")
        self.assertEqual(C(1).is_capacitor, True, "is_capacitor")
        self.assertEqual(L(1).is_inductor, True, "is_inductor")
        self.assertEqual(I(1).is_current_source, True, "is_current_source")
        self.assertEqual(V(1).is_voltage_source, True, "is_voltage_source")

    def test_subs(self):

        n = (R('R1') + C('C1') + L('L1')) | C('C0')
        m = n.subs({'R1':2, 'C1':3, 'L1':4, 'C0':5})
        p = (R(2) + C(3) + L(4)) | C(5)

        self.assertEqual(m.Z(s), p.Z(s), "Z")
        
    def test_noisy(self):

        n = R('R1') + R('R2')        
        m = n.noisy()

        Vn = m.Voc.n
        Vn2 = expr('sqrt(4 * k_B * T * (R1 + R2))')

        n = G('1/R1') + G('1/R2')        
        m = n.noisy()
        Vn = m.Voc.n        
        
        self.assertEqual(Vn, Vn2, "G noise sum")
        
