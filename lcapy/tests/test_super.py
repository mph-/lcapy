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
        self.assertEqual2(Isuper(3) * Zs(2), Vsuper(6), "Isuper(3) * Zs(2)")
        self.assertEqual2(Isuper(12) / Ys(2), Vsuper(6), "Isuper(12) / Ys(2)")

        
    def test_Vsuper_mul_div(self):
        self.assertEqual2(Vsuper(3) * Ys(2), Isuper(6), "Vsuper(3) * Ys(2)")
        self.assertEqual2(Vsuper(12) / Zs(2), Isuper(6), "Vsuper(12) / Zs(2)")

    def test_Vsuper_noise(self):
        self.assertEqual((Vn(3) + Vn(4)).expr, Vn(5).expr, "Vn(3) + Vn(4)")
        self.assertEqual((Vsuper(Vn(3)) + Vsuper(Vn(4))).n.expr,
                          Vsuper(Vn(5)).n.expr,
                          "Vsuper(Vn(3)) + Vsuper(Vn(4))")
        
        
        
        
