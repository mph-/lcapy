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

    def test_Vsuper_dc(self):

        self.assertEqual2(Vsuper(3).dc, 3, "Vsuper(3).dc")
        self.assertEqual2(Vsuper(2, 3).dc, 5, "Vsuper(2, 3).dc")
        self.assertEqual2(Vsuper(2, 3).ac, 0, "Vsuper(2, 3).ac")
        self.assertEqual2(-Vsuper(2).dc, -2, "-Vsuper(2).dc")
        self.assertEqual2(Vsuper(2) + Vsuper(3), Vsuper(5),
                          "Vsuper(2) + Vsuper(3)")
        self.assertEqual2(Vsuper(2) - Vsuper(3), Vsuper(-1),
                          "Vsuper(2) - Vsuper(3)")

    def test_Isuper_dc(self):

        self.assertEqual2(Isuper(3).dc, 3, "Isuper(3).dc")
        self.assertEqual2(Isuper(2, 3).dc, 5, "Isuper(2, 3).dc")
        self.assertEqual2(Isuper(2, 3).ac, 0, "Isuper(2, 3).ac")
        self.assertEqual2(-Isuper(2).dc, -2, "-Isuper(2).dc")
        self.assertEqual2(Isuper(2) + Isuper(3), Isuper(5),
                          "Isuper(2) + Isuper(3)")
        self.assertEqual2(Isuper(2) - Isuper(3), Isuper(-1),
                          "Isuper(2) - Isuper(3)")
        
        
        
        
