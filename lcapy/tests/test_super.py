from lcapy import *
import unittest
import sympy as sym

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

    def test_Voltage_properties(self):
        self.assertEqual(Voltage(3).is_dc, True, "Voltage(3).is_dc")
        self.assertEqual(Voltage(Vphasor(3)).is_ac, True, "Voltage(Vphasor(3)).is_ac")
        self.assertEqual(Voltage(Vconst(2), Vphasor(3)).is_ac, False,
                          "Voltage(Vconst(2), Vphasor(3)).is_ac")
        self.assertEqual(Voltage(Vconst(2), Vphasor(3)).is_ac, False,
                          "Voltage(Vconst(2), Vphasor(3)).is_dc")

    def test_Voltage_add_sub_dc(self):
        self.assertEqual2(Voltage(3).dc, 3, "Voltage(3).dc")
        self.assertEqual2(Voltage(2, 3).dc, 5, "Voltage(2, 3).dc")
        self.assertEqual2(Voltage(2, 3).ac, {}, "Voltage(2, 3).ac")
        self.assertEqual2(-Voltage(2).dc, -2, "-Voltage(2).dc")
        self.assertEqual2(Voltage(2) + Voltage(3), Voltage(5),
                          "Voltage(2) + Voltage(3)")
        self.assertEqual2(Voltage(2) - Voltage(3), Voltage(-1),
                          "Voltage(2) - Voltage(3)")

    def test_Current_add_sub_dc(self):
        self.assertEqual2(Current(3).dc, 3, "Current(3).dc")
        self.assertEqual2(Current(2, 3).dc, 5, "Current(2, 3).dc")
        self.assertEqual2(Current(2, 3).ac, {}, "Current(2, 3).ac")
        self.assertEqual2(-Current(2).dc, -2, "-Current(2).dc")
        self.assertEqual2(Current(2) + Current(3), Current(5),
                          "Current(2) + Current(3)")
        self.assertEqual2(Current(2) - Current(3), Current(-1),
                          "Current(2) - Current(3)")
        
    def test_Current_mul_div(self):
        self.assertEqual2(Current(3) * Impedance(2), Voltage(6), "Current(3) * Impedance(2)")
        self.assertEqual2(Current(12) / Admittance(2), Voltage(6), "Current(12) / Admittance(2)")

        
    def test_Voltage_mul_div(self):
        self.assertEqual2(Voltage(3) * Admittance(2), Current(6), "Voltage(3) * Admittance(2)")
        self.assertEqual2(Voltage(12) / Impedance(2), Current(6), "Voltage(12) / Impedance(2)")

    def test_Voltage_noise(self):
        self.assertEqual((Vn(3) + Vn(4)).expr, Vn(5).expr, "Vn(3) + Vn(4)")
        self.assertEqual((Voltage(Vn(3)) + Voltage(Vn(4))).n.expr,
                          Voltage(Vn(5)).n.expr,
                          "Voltage(Vn(3)) + Voltage(Vn(4))")
        
    def test_Voltage_has(self):

        a = Voltage('3 * exp(-t) * t * a')
        self.assertEqual(a.has(3), True, "has(3)")
        self.assertEqual(a.has(4), False, "has(4)")
        self.assertEqual(a.has(t), True, "has(t)")
        self.assertEqual(a.has_symbol(t), True, "has_symbol(t)")
        self.assertEqual(a.has_symbol('a'), True, "has_symbol(a)")
        self.assertEqual(a.has_symbol('b'), False, "has_symbol(b)")

    def test_Voltage_transform(self):

        V1 = Voltage('3 * exp(-2 * t)')
        self.assertEqual(V1.transform(s), 3 / (s + 2), 'transform(s)')        
        self.assertEqual(V1.transform(jomega), 3 / (j * omega + 2), 'transform(jomega)')

        V2 = Voltage('3 * exp(-2 * t) * u(t)')
        self.assertEqual(V2.transform(s), 3 / (s + 2), 'transform(s)')        
        self.assertEqual(V2.transform(jomega), 3 / (j * omega + 2), 'transform(jomega)')
        self.assertEqual(simplify(V2.transform(f) - 3 / (j * 2 * pi * f + 2)), 0, 'transform(f)')                

        
    def test_Voltage_subs(self):

        a = Voltage('V1')
        b = a.subs('V1', 1)
        c = Voltage(1)
        self.assertEqual(b, c, "Voltage.subs")
