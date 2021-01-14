from lcapy import *
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

    def test_phasor_rms(self):        

        p = phasor(2)
        
        self.assertEqual(p.rms(), sqrt(2), 'rms')
        
        
