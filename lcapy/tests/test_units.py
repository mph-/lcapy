from lcapy import *
from lcapy.units import as_value_unit, u as uu
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_time_units(self):

        self.assertEqual(str(voltage(texpr(1)).units), 'V', 'time voltage')
        self.assertEqual(str(current(texpr(1)).units), 'A', 'time current')
        self.assertEqual(str(admittance(texpr(1)).units), 'S/s', 'time admittance')        
        self.assertEqual(str(impedance(texpr(1)).units), 'ohm/s', 'time impedance')
        self.assertEqual(str((voltage(texpr(1))**2).units), 'V**2', 'time voltage squared')
        self.assertEqual(str((current(texpr(1))**2).units), 'A**2', 'time current squared')
        self.assertEqual(str((admittance(texpr(1))**2).units), 'S**2/s**2', 'time admittance squared')                
        self.assertEqual(str((impedance(texpr(1))**2).units), 'ohm**2/s**2', 'time impedance squared')                

    def test_fourier_units(self):

        self.assertEqual(str(voltage(fexpr(1)).units), 'V/Hz', 'fourier voltage')
        self.assertEqual(str(current(fexpr(1)).units), 'A/Hz', 'fourier current')
        self.assertEqual(str(admittance(fexpr(1)).units), 'S', 'fourier admittance')        
        self.assertEqual(str(impedance(fexpr(1)).units), 'ohm', 'fourier impedance')
        self.assertEqual(str((voltage(fexpr(1))**2).units), 'V**2/Hz**2', 'fourier voltage squared')
        self.assertEqual(str((current(fexpr(1))**2).units), 'A**2/Hz**2', 'fourier current squared')
        self.assertEqual(str((admittance(fexpr(1))**2).units), 'S**2', 'fourier admittance squared')                
        self.assertEqual(str((impedance(fexpr(1))**2).units), 'ohm**2', 'fourier impedance squared')

    def test_constant_units(self):

        self.assertEqual(str(voltage(cexpr(1)).units), 'V', 'constant voltage')
        self.assertEqual(str(current(cexpr(1)).units), 'A', 'constant current')
        self.assertEqual(str(admittance(cexpr(1)).units), 'S', 'constant admittance')        
        self.assertEqual(str(impedance(cexpr(1)).units), 'ohm', 'constant impedance')
        self.assertEqual(str((voltage(cexpr(1))**2).units), 'V**2', 'constant voltage squared')
        self.assertEqual(str((current(cexpr(1))**2).units), 'A**2', 'constant current squared')
        self.assertEqual(str((admittance(cexpr(1))**2).units), 'S**2', 'constant admittance squared')                
        self.assertEqual(str((impedance(cexpr(1))**2).units), 'ohm**2', 'constant impedance squared')

        self.assertEqual(str(resistance(1).units), 'ohm', 'resistance')
        self.assertEqual(str(conductance(1).units), 'S', 'conductance')
        self.assertEqual(str(capacitance(1).units), 'F', 'capacitance')
        self.assertEqual(str(inductance(1).units), 'H', 'inductance')
        self.assertEqual(str(reactance(1).units), 'ohm', 'reactance')
        self.assertEqual(str(susceptance(1).units), 'S', 'susceptance')                
        
    def test_phasor_units(self):

        self.assertEqual(str(voltage(phasor(1)).units), 'V', 'phasor voltage')
        self.assertEqual(str(current(phasor(1)).units), 'A', 'phasor current')
        self.assertEqual(str(admittance(phasor(1)).units), 'S', 'phasor admittance')        
        self.assertEqual(str(impedance(phasor(1)).units), 'ohm', 'phasor impedance')
        self.assertEqual(str((voltage(phasor(1))**2).units), 'V**2', 'phasor voltage squared')
        self.assertEqual(str((current(phasor(1))**2).units), 'A**2', 'phasor current squared')
        self.assertEqual(str((admittance(phasor(1))**2).units), 'S**2', 'phasor admittance squared')                
        self.assertEqual(str((impedance(phasor(1))**2).units), 'ohm**2', 'phasor impedance squared')

    def test_discrete_time_units(self):

        self.assertEqual(str(voltage(nexpr(1)).units), 'V', 'nexpr voltage')
        self.assertEqual(str(current(nexpr(1)).units), 'A', 'nexpr current')
        self.assertEqual(str(admittance(nexpr(1)).units), 'S', 'nexpr admittance')        
        self.assertEqual(str(impedance(nexpr(1)).units), 'ohm', 'nexpr impedance')
        self.assertEqual(str((voltage(nexpr(1))**2).units), 'V**2', 'nexpr voltage squared')
        self.assertEqual(str((current(nexpr(1))**2).units), 'A**2', 'nexpr current squared')
        self.assertEqual(str((admittance(nexpr(1))**2).units), 'S**2', 'nexpr admittance squared')                
        self.assertEqual(str((impedance(nexpr(1))**2).units), 'ohm**2', 'nexpr impedance squared')

    def test_discrete_fourier_units(self):

        self.assertEqual(str(voltage(kexpr(1)).units), 'V', 'kexpr voltage')
        self.assertEqual(str(current(kexpr(1)).units), 'A', 'kexpr current')
        self.assertEqual(str(admittance(kexpr(1)).units), 'S', 'kexpr admittance')        
        self.assertEqual(str(impedance(kexpr(1)).units), 'ohm', 'kexpr impedance')
        self.assertEqual(str((voltage(kexpr(1))**2).units), 'V**2', 'kexpr voltage squared')
        self.assertEqual(str((current(kexpr(1))**2).units), 'A**2', 'kexpr current squared')
        self.assertEqual(str((admittance(kexpr(1))**2).units), 'S**2', 'kexpr admittance squared')
        self.assertEqual(str((impedance(kexpr(1))**2).units), 'ohm**2', 'kexpr impedance squared')
        
    def test_new_units(self):

        v = voltage(3)
        i = current(4)

        p = v * i
        # Convert V * W to W if canonical_units not enabled.
        p = p.simplify_units()

        self.assertEqual(str(p.units), 'W', 'power')        

    def test_functions(self):

        self.assertEqual(str(delta(t).units), '1/s', 'delta(t)')
        self.assertEqual(str(diff(t, t).units), '1', 'diff(t, t)')
        self.assertEqual(str(integrate(t, t).units), 's**2', 'integrate(t, t)')
        self.assertEqual(str(integrate(t, (t, 0, 1)).units), 's**2', 'integrate(t, (t, 0, 1))')
        self.assertEqual(str(atan(f).units), 'rad', 'atan(f)')        
        
    def test_convolution(self):

        i = current('i(t)')
        z = impedance('z(t)')
        v = i.convolve(z).simplify_units()
        
        self.assertEqual(str(v.units), 'V', 'i convolve z')                

    def test_as_value_unit(self):

        i = current(7)
        val, unit = as_value_unit(i.expr_with_units)
        self.assertEqual(str(unit), 'A', 'test unit')
        self.assertEqual(val, 7, 'test val')
