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

    def test_constant_units(self):

        self.assertEqual(voltage(cexpr(1)).units, 'V', 'constant voltage')
        self.assertEqual(current(cexpr(1)).units, 'A', 'constant current')
        self.assertEqual(admittance(cexpr(1)).units, 'S', 'constant admittance')        
        self.assertEqual(impedance(cexpr(1)).units, 'ohm', 'constant impedance')
        self.assertEqual((voltage(cexpr(1))**2).units, 'V^2', 'constant voltage squared')
        self.assertEqual((current(cexpr(1))**2).units, 'A^2', 'constant current squared')
        self.assertEqual((admittance(cexpr(1))**2).units, 'S^2', 'constant admittance squared')                
        self.assertEqual((impedance(cexpr(1))**2).units, 'ohm^2', 'constant impedance squared')                
        
    def test_phasor_units(self):

        self.assertEqual(voltage(phasor(1)).units, 'V', 'phasor voltage')
        self.assertEqual(current(phasor(1)).units, 'A', 'phasor current')
        self.assertEqual(admittance(phasor(1)).units, 'S', 'phasor admittance')        
        self.assertEqual(impedance(phasor(1)).units, 'ohm', 'phasor impedance')
        self.assertEqual((voltage(phasor(1))**2).units, 'V^2', 'phasor voltage squared')
        self.assertEqual((current(phasor(1))**2).units, 'A^2', 'phasor current squared')
        self.assertEqual((admittance(phasor(1))**2).units, 'S^2', 'phasor admittance squared')                
        self.assertEqual((impedance(phasor(1))**2).units, 'ohm^2', 'phasor impedance squared')

    def test_discrete_time_units(self):

        self.assertEqual(voltage(nexpr(1)).units, 'V', 'nexpr voltage')
        self.assertEqual(current(nexpr(1)).units, 'A', 'nexpr current')
        self.assertEqual(admittance(nexpr(1)).units, 'S', 'nexpr admittance')        
        self.assertEqual(impedance(nexpr(1)).units, 'ohm', 'nexpr impedance')
        self.assertEqual((voltage(nexpr(1))**2).units, 'V^2', 'nexpr voltage squared')
        self.assertEqual((current(nexpr(1))**2).units, 'A^2', 'nexpr current squared')
        self.assertEqual((admittance(nexpr(1))**2).units, 'S^2', 'nexpr admittance squared')                
        self.assertEqual((impedance(nexpr(1))**2).units, 'ohm^2', 'nexpr impedance squared')

    def test_discrete_fourier_units(self):

        self.assertEqual(voltage(kexpr(1)).units, 'V', 'kexpr voltage')
        self.assertEqual(current(kexpr(1)).units, 'A', 'kexpr current')
        self.assertEqual(admittance(kexpr(1)).units, 'S', 'kexpr admittance')        
        self.assertEqual(impedance(kexpr(1)).units, 'ohm', 'kexpr impedance')
        self.assertEqual((voltage(kexpr(1))**2).units, 'V^2', 'kexpr voltage squared')
        self.assertEqual((current(kexpr(1))**2).units, 'A^2', 'kexpr current squared')
        self.assertEqual((admittance(kexpr(1))**2).units, 'S^2', 'kexpr admittance squared')                
        self.assertEqual((impedance(kexpr(1))**2).units, 'ohm^2', 'kexpr impedance squared')                        
