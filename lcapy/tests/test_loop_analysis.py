from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_loop(self):

        a = Circuit("""
        V 1 0 {v(t)}
        R1 1 2
        R2 2 3
        R3 3 0_3
        W 0 0_3
        W 3 3_a
        R4 3_a 0_4
        W 0_3 0_4""")

        la = LoopAnalysis(a)

        la_eqs = la.mesh_equations()
        loops = la.loops_by_cpt_name()

        key = list(la_eqs.keys())[0]        
        eq = la_eqs[key]
        loop = loops[0]
        
        if loop == ['R2', 'C']:
            self.assertEqual(eq.lhs, voltage('R3*(I_2(t) - I_1(t)) + R4*(I_2(t) - I_1(t)) + v(t)'), 'mesh_equations()[1].lhs')
        elif loop == ['V', 'R1', 'R2', 'R3']:
            self.assertEqual(eq.lhs, voltage('-R1*I_1(t) - R2*I_1(t) + R3*(-I_1(t) + I_2(t)) + v(t)'), 'mesh_equations()[1].lhs')
        elif loop == ['R3', 'R4']:
            self.assertEqual(eq.lhs, voltage('R3*(I_1(t) - I_2(t)) + R4*(I_1(t) - I_2(t))'), 'mesh_equations()[1].lhs')            
        else:
            raise ValueError('Unhandled loop %s' % loop)

        # Use abs since sometimes the current flows in the other direction
        # depending how the loop was chosen.

        self.assertEqual(eq.rhs, voltage(0), 'mesh_equations()[1].rhs')
        
        
