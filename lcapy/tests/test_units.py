from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_time_units(self):

        self.assertEqual(voltage(texpr(1)).units, 'V', 'time voltage')
        self.assertEqual(current(texpr(1)).units, 'A', 'time current')
        self.assertEqual(admittance(texpr(1)).units, 'S/s', 'time admittance')        
        self.assertEqual(impedance(texpr(1)).units, 'ohm/s', 'time impedance')
        self.assertEqual((voltage(texpr(1))**2).units, 'V^2', 'time voltage squared')
        self.assertEqual((current(texpr(1))**2).units, 'A^2', 'time current squared')
        self.assertEqual((admittance(texpr(1))**2).units, 'S^2/s^2', 'time admittance squared')                
        self.assertEqual((impedance(texpr(1))**2).units, 'ohm^2/s^2', 'time impedance squared')                

    def test_fourier_units(self):

        self.assertEqual(voltage(fexpr(1)).units, 'V/Hz', 'fourier voltage')
        self.assertEqual(current(fexpr(1)).units, 'A/Hz', 'fourier current')
        self.assertEqual(admittance(fexpr(1)).units, 'S', 'fourier admittance')        
        self.assertEqual(impedance(fexpr(1)).units, 'ohm', 'fourier impedance')
        self.assertEqual((voltage(fexpr(1))**2).units, 'V^2/Hz^2', 'fourier voltage squared')
        self.assertEqual((current(fexpr(1))**2).units, 'A^2/Hz^2', 'fourier current squared')
        self.assertEqual((admittance(fexpr(1))**2).units, 'S^2', 'fourier admittance squared')                
        self.assertEqual((impedance(fexpr(1))**2).units, 'ohm^2', 'fourier impedance squared')        
        
