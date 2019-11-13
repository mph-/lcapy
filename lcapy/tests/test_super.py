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

    def test_Vsuper_properties(self):
        self.assertEqual(Vsuper(3).is_dc, True, "Vsuper(3).is_dc")
        self.assertEqual(Vsuper(Vphasor(3)).is_ac, True, "Vsuper(Vphasor(3)).is_ac")
        self.assertEqual(Vsuper(Vconst(2), Vphasor(3)).is_ac, False,
                          "Vsuper(Vconst(2), Vphasor(3)).is_ac")
        self.assertEqual(Vsuper(Vconst(2), Vphasor(3)).is_ac, False,
                          "Vsuper(Vconst(2), Vphasor(3)).is_dc")

    def test_Vsuper_add_sub_dc(self):
        self.assertEqual2(Vsuper(3).dc, 3, "Vsuper(3).dc")
        self.assertEqual2(Vsuper(2, 3).dc, 5, "Vsuper(2, 3).dc")
        self.assertEqual2(Vsuper(2, 3).ac, {}, "Vsuper(2, 3).ac")
        self.assertEqual2(-Vsuper(2).dc, -2, "-Vsuper(2).dc")
        self.assertEqual2(Vsuper(2) + Vsuper(3), Vsuper(5),
                          "Vsuper(2) + Vsuper(3)")
        self.assertEqual2(Vsuper(2) - Vsuper(3), Vsuper(-1),
                          "Vsuper(2) - Vsuper(3)")

    def test_Isuper_add_sub_dc(self):
        self.assertEqual2(Isuper(3).dc, 3, "Isuper(3).dc")
        self.assertEqual2(Isuper(2, 3).dc, 5, "Isuper(2, 3).dc")
        self.assertEqual2(Isuper(2, 3).ac, {}, "Isuper(2, 3).ac")
        self.assertEqual2(-Isuper(2).dc, -2, "-Isuper(2).dc")
        self.assertEqual2(Isuper(2) + Isuper(3), Isuper(5),
                          "Isuper(2) + Isuper(3)")
        self.assertEqual2(Isuper(2) - Isuper(3), Isuper(-1),
                          "Isuper(2) - Isuper(3)")
        
    def test_Isuper_mul_div(self):
        self.assertEqual2(Isuper(3) * ZZ(2), Vsuper(6), "Isuper(3) * ZZ(2)")
        self.assertEqual2(Isuper(12) / YY(2), Vsuper(6), "Isuper(12) / YY(2)")

        
    def test_Vsuper_mul_div(self):
        self.assertEqual2(Vsuper(3) * YY(2), Isuper(6), "Vsuper(3) * YY(2)")
        self.assertEqual2(Vsuper(12) / ZZ(2), Isuper(6), "Vsuper(12) / ZZ(2)")

    def test_Vsuper_noise(self):
        self.assertEqual((Vn(3) + Vn(4)).expr, Vn(5).expr, "Vn(3) + Vn(4)")
        self.assertEqual((Vsuper(Vn(3)) + Vsuper(Vn(4))).n.expr,
                          Vsuper(Vn(5)).n.expr,
                          "Vsuper(Vn(3)) + Vsuper(Vn(4))")
        
    def test_Vsuper_has(self):

        a = Vsuper('3 * exp(-t) * t * a')
        self.assertEqual(a.has(3), True, "has(3)")
        self.assertEqual(a.has(4), False, "has(4)")
        self.assertEqual(a.has(t), True, "has(t)")
        self.assertEqual(a.has_symbol(t), True, "has_symbol(t)")
        self.assertEqual(a.has_symbol('a'), True, "has_symbol(a)")
        self.assertEqual(a.has_symbol('b'), False, "has_symbol(b)")

    def test_Vsuper_transform(self):

        V1 = Vsuper('3 * exp(-2 * t)')
        self.assertEqual(V1.transform(s), 3 / (s + 2), 'transform(s)')        
        self.assertEqual(V1.transform(jomega), 3 / (j * omega + 2), 'transform(jomega)')

        V2 = Vsuper('3 * exp(-2 * t) * u(t)')
        self.assertEqual(V2.transform(s), 3 / (s + 2), 'transform(s)')        
        self.assertEqual(V2.transform(jomega), 3 / (j * omega + 2), 'transform(jomega)')
        self.assertEqual(simplify(V2.transform(f) - 3 / (j * 2 * pi * f + 2)), 0, 'transform(f)')                

        
    def test_Vsuper_subs(self):

        a = Vsuper('V1')
        b = a.subs('V1', 1)
        c = Vsuper(1)
        self.assertEqual(b, c, "Vsuper.subs")
