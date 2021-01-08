from lcapy import *
from lcapy.discretetime import *
import unittest


class LcapyTesterLatex(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_latex_subscripts(self):

        self.assertEqual(expr('V1').latex(), 'V_{1}', 'V1')
        self.assertEqual(expr('V_1').latex(), 'V_{1}', 'V_1')
        self.assertEqual(expr('V_C').latex(), 'V_{C}', 'V_C')
        
        self.assertEqual(expr('V_C1').latex(), 'V_{C_{1}}', 'V_C1')
        self.assertEqual(expr('V_C_1').latex(), 'V_{C_{1}}', 'V_C_1')
        
        self.assertEqual(expr('Vrms').latex(), 'V_{\\mathrm{rms}}', 'Vrms')
        self.assertEqual(expr('V_rms').latex(), 'V_{\\mathrm{rms}}', 'V_rms')

        self.assertEqual(expr('V_C1(t)').latex(), 'V_{C_{1}}(t)', 'V_C1(t)')
        self.assertEqual(expr('V_C_1(t)').latex(), 'V_{C_{1}}(t)', 'V_C_1(t)')        

    def test_latex_greek(self):

        self.assertEqual(expr('alpha').latex(), '\\alpha', 'alpha')        

