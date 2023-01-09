from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_loop1(self):

        a = Circuit("""
        V 1 0 {v(t)}
        R 1 0""")

        la = a.loop_analysis()

        la_eqs = la.mesh_equations()
        loops = la.loops_by_cpt_name()

        key = list(la_eqs.keys())[0]
        eq = la_eqs[key]

        self.assertEqual(len(loops), 1, 'num loops')

        loop = loops[0]

        self.assertEqual(eq.lhs, voltage('v(t) - R * i_1(t)'),
                         'mesh_equations()[1].lhs')

        self.assertEqual(eq.rhs, voltage(0), 'mesh_equations()[1].rhs')

    def test_loop2(self):

        a = Circuit("""
        V1 1 0 {v1(t)}
        V2 1 2 {v2(t)}
        R 2 0""")

        la = a.loop_analysis()

        la_eqs = la.mesh_equations()
        loops = la.loops_by_cpt_name()

        key = list(la_eqs.keys())[0]
        eq = la_eqs[key]

        self.assertEqual(len(loops), 1, 'num loops')

        loop = loops[0]

        self.assertEqual(eq.lhs,
                         voltage('v1(t) - v2(t) - R * i_1(t)'),
                         'mesh_equations()[1].lhs')

        self.assertEqual(eq.rhs, voltage(0), 'mesh_equations()[1].rhs')

    def test_loop3(self):

        a = Circuit("""
        V 1 0 {v(t)}
        R1 1 2
        R2 2 3
        R3 3 0_3
        W 0 0_3
        W 3 3_a
        R4 3_a 0_4
        W 0_3 0_4""")

        la = a.loop_analysis()

        la_eqs = la.mesh_equations()
        loops = la.loops_by_cpt_name()

        key = list(la_eqs.keys())[0]
        eq = la_eqs[key]
        loop = loops[0]

        if loop == ['R2', 'C']:
            self.assertEqual(abs(eq.lhs), voltage(
                'abs(R3*(i_2(t) - i_1(t)) + R4*(i_2(t) - i_1(t)) + v(t))'), 'mesh_equations()[1].lhs')
        elif loop == ['V', 'R1', 'R2', 'R3']:
            self.assertEqual(abs(eq.lhs), voltage(
                'abs(-R1*i_1(t) - R2*i_1(t) + R3*(-i_1(t) + i_2(t)) + v(t))'), 'mesh_equations()[1].lhs')
        elif loop == ['R3', 'R4']:
            self.assertEqual(abs(eq.lhs), voltage(
                'abs(R3*(i_1(t) - i_2(t)) + R4*(i_1(t) - i_2(t)))'), 'mesh_equations()[1].lhs')
        else:
            raise ValueError('Unhandled loop %s' % loop)

        # Use abs since sometimes the current flows in the other direction
        # depending how the loop was chosen.

        self.assertEqual(eq.rhs, voltage(0), 'mesh_equations()[1].rhs')
