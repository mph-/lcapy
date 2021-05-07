from lcapy import *
from lcapy.phasor import PhasorFrequencyDomainExpression
import unittest

class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy"""

    def assertEqual2(self, ans1, ans2, comment):

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            pprint(ans1)
            pprint(ans2)
            raise AssertionError(e)

    def test_phasor(self):

        a = phasor(3, omega=7)
        self.assertEqual(a.omega, 7, 'omega')
        self.assertEqual(a.phase, 0, 'phase')
        self.assertEqual(a.magnitude, 3, 'magnitude')
        self.assertEqual(a.N, 3, 'N')
        self.assertEqual(a.D, 1, 'D')                        
        
        a = phasor(-3, omega=7)
        self.assertEqual(a.phase, pi, 'phase')                
        self.assertEqual(a.magnitude, 3, 'magnitude')

        a = phasor(-3 + 4j, omega=7)
        self.assertEqual(a.magnitude, 5, 'magnitude')                        

        v = voltage('cos(5 * t)')
        i = current('cos(5 * t)')
        V = v.as_phasor()
        I = i.as_phasor()

        Z = V / I
        self.assertEqual(Z.is_impedance, True, 'impedance')
        self.assertEqual(Z.omega, 5, 'omega')

        V2 = phasorvoltage('cos(5 * t)')
        I2 = phasorcurrent('cos(5 * t)')
        self.assertEqual(V, V2, 'phasorvoltage')
        self.assertEqual(I, I2, 'phasorcurrent')                
        
    def test_phasor_rms(self):        

        p = phasor(2)
        
        self.assertEqual(p.rms(), sqrt(2), 'rms')

    def test_phasor_mul_omega(self):

        p = phasor(2)
        q = p * omega
        r = phasor(2 * omega)

        self.assertEqual(q, r, 'phasor(2) * omega')        
        
    def test_phasor_transform(self):

        H = 1 / s
        P = H(jw)

        self.assertEqual(type(P), PhasorFrequencyDomainExpression, 'H(jw)')

        P = H(j * 2 * pi * f)
        self.assertEqual(type(P), PhasorFrequencyDomainExpression, 'H(j * 2 * pi * f)')

        self.assertEqual(P.var, None, 'P.var')        
        
        
        
