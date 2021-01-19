from lcapy import *
from lcapy.texpr import TimeDomainExpression
from lcapy.noiseomegaexpr import AngularFourierNoiseDomainVoltage
import unittest


class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy noise"""

    def assertEqual2(self, ans1, ans2, comment):

        ans1 = ans1.canonical()
        ans2 = ans2.canonical()

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            ans1.pprint()
            ans2.pprint()
            raise AssertionError(e)

    def test_noise(self):
        a = AngularFourierNoiseDomainVoltage(2)
        self.assertEqual(a.nid, a.conjugate.nid, "Different nids for conjugate")
        self.assertEqual(a.nid, a.real.nid, "Different nids for real")
        self.assertEqual(a.nid, a.imag.nid, "Different nids for imag")
        
    def test_noise1(self):
        """Lcapy: check circuit noise for voltage divider"""

        a = Circuit()
        a.add('V1 1 0 noise 3') 
        a.add('R1 1 2 2')
        a.add('R2 2 0 4')
        V1 = a.R1.V.n
        self.assertEqual2(V1, AngularFourierNoiseDomainVoltage(1, nid=V1.nid), "Incorrect ratio")

    def test_noise2(self):
        """Lcapy: check circuit noise for pair of sources"""

        a = Circuit()
        a.add('V1 1 0 noise 3')
        a.add('V2 2 1 noise 4')
        a.add('R1 2 0 5')
        V1 = a.R1.V.n
        self.assertEqual2(V1, AngularFourierNoiseDomainVoltage(5, nid=V1.nid), "Incorrect noise sum")        
        
    def test_filtered_noise1(self):
        """Lcapy: check circuit filtered noise"""

        a = Circuit()
        a.add('V1 1 0 noise 3') 
        a.add('R1 1 2 2')
        a.add('C1 2 0 4')         
#        self.assertEqual2(a.R1.V.n, Vnoisy(1), "Incorrect ratio")


    def test_filtered_noise2(self):
        """Lcapy: check circuit filtered noise"""

        a = Circuit()
        a.add('V1 1 0 noise {sqrt(4 * k_B * T * R)}') 
        a.add('R1 1 2 R')
        a.add('C1 2 0 C')
        self.assertEqual2(a.C1.V.n.rms(), voltage('sqrt(k_B * T / C)'),
                          "Incorrect capacitor voltage")        

    def test_filtered_noise3(self):
        """Lcapy: check circuit filtered noise"""

        a = Circuit()
        a.add('V1 1 0 noise 20') 
        a.add('R1 1 2 1')
        a.add('C1 2 0 2')         
        self.assertEqual(a.C1.V.n.rms(), voltage(5 * sqrt(2)),
                         "Incorrect capacitor voltage")

    def test_noisy1(self):

        a = Circuit()
        a.add('R1 1 0')
        a.add('R2 1 0')
        an = a.noisy()
        b = Circuit()
        b.add('R1 1 0 {R1 * R2 / (R1 + R2)}')
        bn = b.noisy()
        self.assertEqual(an[1].V.n.expr, bn[1].V.n.expr, "Incorrect noise")

    def test_noisy_transform1(self):

        a = Circuit()
        a.add('R1 1 0')
        a.add('R2 1 0')
        an = a.noisy()
        Vn1 = an.R1.V.n

        Vn2 = Vn1(f)(omega)
        
        self.assertEqual(Vn1.nid, Vn2.nid, "Incorrect noise nid")
        self.assertEqual(Vn1, Vn2, "Incorrect noise nid")

    def test_noisy_transform_add1(self):

        a = AngularFourierNoiseDomainVoltage(2 * omega)
        b = AngularFourierNoiseDomainVoltage(3 + omega)
        
    def test_noisef_to_noiseomega(self):

        a = AngularFourierNoiseDomainVoltage(omega)
        b = a(omega)
        c = a(f)
        d = c(omega)

        self.assertEqual(a, b, 'noisyomega -> noisyomega')
        self.assertEqual(a, d, 'noisyomega -> noisyf -> noisyomega')
        
    def test_noisevoltage(self):

        V1 = noisevoltage(1)
        self.assertEqual(V1.is_fourier_noise_domain, True, 'check domain')
        self.assertEqual(V1.is_voltage, True, 'check quantity')

        V2 = noisevoltage(1 / f)
        self.assertEqual(V2.is_fourier_noise_domain, True, 'check domain')

        V3 = noisevoltage(1 / omega)
        self.assertEqual(V3.is_angular_fourier_noise_domain, True, 'check domain')        

        x = V1.sample(4)
        self.assertEqual(len(x), 4, 'sample size')                
        
    def test_noisecurrent(self):

        I1 = noisecurrent(1)
        self.assertEqual(I1.is_fourier_noise_domain, True, 'check domain')
        self.assertEqual(I1.is_current, True, 'check quantity')

        I2 = noisecurrent(1 / f)
        self.assertEqual(I2.is_fourier_noise_domain, True, 'check domain')

        I3 = noisecurrent(1 / omega)
        self.assertEqual(I3.is_angular_fourier_noise_domain, True, 'check domain')        

    def test_noisevoltage_init(self):

        V1 = noisevoltage(1)
        # Check that nid copied
        V2 = V1.as_voltage()

        self.assertEqual(V1, V2, 'as_voltage')        
