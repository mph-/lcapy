from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_nodal(self):

        a = Circuit("""
        V 1 0 {v(t)}
        R1 1 2
        L 2 3
        R2 3 0_3
        W 0 0_3
        W 3 3_a
        C 3_a 0_4
        W 0_3 0_4""")

        na = NodalAnalysis(a)

        na_eqs = na.nodal_equations()

        self.assertEqual(na_eqs[1].lhs, expr('V1(t)'), 'nodal_equations()[1].lhs')
        self.assertEqual(na_eqs[1].rhs, expr('v(t)'), 'nodal_equations()[1].rhs')        
        
