from lcapy import *
from lcapy.phasor import PhasorDomainExpression
from lcapy.sexpr import LaplaceDomainVoltage
import unittest
import sympy as sym
from lcapy.sym import omega0sym

class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def assertEqual2(self, ans1, ans2, comment):

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            ans1.pprint()
            ans2.pprint()
            raise AssertionError(e)

    def test_R(self):
        """Lcapy: check R

        """
        a = R(10)
        self.assertEqual(a.Z, 10, "Z incorrect.")
        self.assertEqual(a.Y * 10, 1, "Y incorrect.")        
        self.assertEqual2(a.Voc.s, 0, "Voc incorrect.")
        self.assertEqual2(a.Isc.s, 0, "Isc incorrect.")
        self.assertEqual2(a.voc, 0, "voc incorrect.")
        self.assertEqual2(a.isc, 0, "isc incorrect.")

    def test_L(self):
        """Lcapy: check L

        """
        a = L(10, 5)
        self.assertEqual2(a.Z(s), 10 * s, "Zs incorrect.")
        self.assertEqual2(a.Y(s), 1 / (10 * s), "Ys incorrect.")
        self.assertEqual2(a.Z(jw), 10 * j * omega, "Z incorrect.")
        self.assertEqual2(a.Y(jw), 1 / (10 * j * omega), "Y incorrect.")                
        self.assertEqual2(a.Voc.s, -10 * 5, "Voc incorrect.")
        self.assertEqual2(a.Isc.s, -5 / s, "Isc incorrect.")        

    def test_C(self):
        """Lcapy: check C

        """
        a = C(10, 5)

        self.assertEqual2(a.Z(s), 1 / (10 * s), "Zs incorrect.")
        self.assertEqual2(a.Y(s), 10 * s, "Ys incorrect.")
        self.assertEqual2(a.Z(jw), 1 / (10 * j * omega), "Z incorrect.")
        self.assertEqual2(a.Y(jw), 10 * j * omega, "Y incorrect.")                
        self.assertEqual2(a.Voc.s, 5 / s, "Voc incorrect.")
        self.assertEqual(a.Isc.s, 50, "Isc incorrect.")        

    def test_R_series_R(self):
        """Lcapy: check R + R

        """
        a = R(10) + R(5)
        b = a.simplify()
        self.assertEqual2(b.R, 15, "R incorrect.")
        self.assertEqual2(b.Z, 15, "Z incorrect.")
        self.assertEqual(b.Voc.s, 0, "Voc incorrect.")
        self.assertEqual(b.Isc.s, 0, "Isc incorrect.")                
        self.assertEqual2(type(b), R, "type incorrect.")

    def test_L_series_L(self):
        """Lcapy: check L + L

        """
        a = L(10) + L(5)
        b = a.simplify()
        self.assertEqual2(b.L, 15, "L incorrect.")
        self.assertEqual2(b.Z(s), 15 * s, "Zs incorrect.")
        self.assertEqual2(b.Z(jw), 15 * j * omega, "Z incorrect.")        
        self.assertEqual2(type(b), L, "type incorrect.")

    def test_C_series_C(self):
        """Lcapy: check C + C

        """
        a = C(10) + C(15)
        b = a.simplify()
        self.assertEqual2(b.C, 6, "C incorrect.")
        self.assertEqual2(b.Z(s), 1 / (6 * s), "Zs incorrect.")
        self.assertEqual2(b.Z(jw), 1 / (6 * j * omega), "Z incorrect.")        
        self.assertEqual2(type(b), C, "type incorrect.")

    def test_Vdc_series_Vdc(self):
        """Lcapy: check Vdc + Vdc

        """
        a = Vdc(10) + Vdc(5)
        b = a.simplify()
        self.assertEqual2(b.voc, Vdc(15).voc, "voc incorrect.")
        self.assertEqual(b.Voc.dc, 15, "Voc incorrect.")
        self.assertEqual2(type(b), Vdc, "type incorrect.")
        self.assertEqual2(b.voc, 15, "voc incorrect.")

    def test_Vdc_series_R(self):
        """Lcapy: check Vdc + R

        """
        a = Vdc(10) + R(5)
        self.assertEqual2(a.voc, Vdc(10).voc, "voc incorrect.")
        self.assertEqual(a.Voc.dc, 10, "Voc incorrect.")
        self.assertEqual(a.V.dc, 10, "Voc incorrect.")
        self.assertEqual(a.I, 0, "I incorrect.")
        self.assertEqual(a.v(0), 10, "Voc incorrect.")                
        self.assertEqual(a.i(0), 0, "I incorrect.")                

        b = a.norton()
        self.assertEqual2(b.voc, Vdc(10).voc, "voc incorrect.")
        self.assertEqual(b.Voc.dc, 10, "Voc incorrect.")
        self.assertEqual(b.V.dc, 10, "Voc incorrect.")
        self.assertEqual(b.I, 0, "I incorrect.")
        self.assertEqual(b.v(0), 10, "Voc incorrect.")                
        self.assertEqual(b.i(0), 0, "I incorrect.")

    def test_Vac_series_R(self):
        """Lcapy: check Vac + R

        """
        a = Vac(10) + R(5)
        self.assertEqual2(a.voc, Vac(10).voc, "voc incorrect.")
        self.assertEqual(a.Voc.ac[omega0sym], PhasorDomainExpression(10), "Voc incorrect.")
        self.assertEqual(a.V.ac[omega0sym], PhasorDomainExpression(10), "Voc incorrect.")
        self.assertEqual(a.I, 0, "I incorrect.")
        self.assertEqual(a.v(0), 10, "Voc incorrect.")                
        self.assertEqual(a.i(0), 0, "I incorrect.")                

        b = a.norton()
        self.assertEqual2(b.voc, Vac(10).voc, "voc incorrect.")
        self.assertEqual(b.Voc.ac[omega0sym], PhasorDomainExpression(10), "Voc incorrect.")
        self.assertEqual(b.V.ac[omega0sym], PhasorDomainExpression(10), "Voc incorrect.")
        self.assertEqual(b.I, 0, "I incorrect.")
        self.assertEqual(b.v(0), 10, "Voc incorrect.")                
        self.assertEqual(b.i(0), 0, "I incorrect.")        
        
    def test_R_series_L(self):
        """Lcapy: check R + L

        """
        a = R(10) + L(5)
        self.assertEqual2(a.Z(s), 5 * s + 10, "Zs incorrect.")

    def test_R_parallel_R(self):
        """Lcapy: check R | R

        """
        a = R(10) | R(15)
        b = a.simplify()
        self.assertEqual2(b.R, 6, "R incorrect.")
        self.assertEqual2(b.Z, 6, "Z incorrect.")
        self.assertEqual2(type(b), R, "type incorrect.")

    def test_L_parallel_L(self):
        """Lcapy: check L | L

        """
        a = L(10) | L(15)
        b = a.simplify()
        self.assertEqual2(b.L, 6, "L incorrect.")
        self.assertEqual2(b.Z(s), 6 * s, "Zs incorrect.")
        self.assertEqual2(type(b), L, "type incorrect.")

    def test_C_parallel_C(self):
        """Lcapy: check C | C

        """
        a = C(10) | C(15)
        b = a.simplify()
        self.assertEqual2(b.C, 25, "C incorrect.")
        self.assertEqual2(b.Z(s), 1 / (25 * s), "Zs incorrect.")
        self.assertEqual2(type(b), C, "type incorrect.")

    def test_I_parallel_I(self):
        """Lcapy: check I | I

        """
        a = Idc(10) | Idc(5)
        b = a.simplify()
        self.assertEqual2(b.isc, Idc(15).isc, "isc incorrect.")
        self.assertEqual2(b.Isc.dc, 15, "Isc incorrect.")
        self.assertEqual2(type(b), Idc, "type incorrect.")
        self.assertEqual2(b.isc, 15, "isc incorrect.")

    def test_load(self):
        """Lcapy: check load

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.load(R(30))

        self.assertEqual2(type(b), Ser, "type incorrect.")
        self.assertEqual2(b.Z, 7.5, "Shunt loaded R incorrect Z.")

    def test_open_circuit(self):
        """Lcapy: check open_circuit

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.open_circuit()

        self.assertEqual2(b.Z, 10, "incorrect Z.")
        self.assertEqual2(b.Voc.laplace(), LaplaceDomainVoltage(5) / s, "incorrect V.")

    def test_short_circuit(self):
        """Lcapy: check short_circuit

        """

        a = Series(R(10) + Vdc(5))
        b = a.short_circuit()

        self.assertEqual2(b.Z, 10, "incorrect Z.")
        self.assertEqual2(b.Voc.laplace(), LaplaceDomainVoltage(5) / s, "incorrect V.")

    def test_v(self):
        """Lcapy: check inverse Laplace for voltage sources"""

        omega0 = symbol('omega_0', real=True)

        a = Vdc(10)
        self.assertEqual(a.voc, 10, "DC incorrect.")        

        a = Vstep(10)
        self.assertEqual2(a.voc, 10 * Heaviside(t), "Step incorrect.")

        a = Vac(10)
        self.assertEqual2(a.voc, 10 * cos(omega0 * t), "AC incorrect.")


    def test_i(self):
        """Lcapy: check inverse Laplace for current sources"""

        omega0 = symbol('omega_0', real=True)

        a = Idc(10)
        self.assertEqual(a.isc, 10, "DC incorrect.")        

        a = Istep(10)
        self.assertEqual2(a.isc, 10 * Heaviside(t), "Step incorrect.")

        a = Iac(10)
        self.assertEqual2(a.isc, 10 * cos(omega0 * t), "AC incorrect.")


    def test_CPE(self):
        """Lcapy: check CPE

        """
        a = CPE(10, 1)

        self.assertEqual2(a.Z(s), 1 / (10 * s), "Zs incorrect.")
        self.assertEqual2(a.Y(s), 10 * s, "Ys incorrect.")
        self.assertEqual2(a.Z(jw), 1 / (10 * j * omega), "Zw incorrect.")
        self.assertEqual2(a.Y(jw), 10 * j * omega, "Yw incorrect.")

        a = CPE(10, 0)

        self.assertEqual2(a.Z(s), 1.0 / 10, "Zs incorrect.")
        self.assertEqual2(a.Y(s), 10, "Ys incorrect.")
        self.assertEqual2(a.Z(jw), 1.0 / 10, "Zw incorrect.")
        self.assertEqual2(a.Y(jw), 10, "Yw incorrect.")

        a = CPE(10, 0.5)

        self.assertEqual2(a.Y(s), 10 * sqrt(s), "Ys incorrect.")

    def simplify1(self):

        a = I(0) + Y(1)
        b = a.simplify()

        self.assertEqual2(b._Z, None, "Z incorrect.")
        self.assertEqual2(b._Y, 1, "Y incorrect.")                
