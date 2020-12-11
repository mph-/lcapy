from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_fourier(self):

        self.assertEqual(DiracDelta(t).fourier(), 1, "DiracDelta(t)")
        self.assertEqual(TimeDomainVoltage('x(t)').fourier(), FourierDomainVoltage('X(f)'), "x(t)")
        self.assertEqual(TimeDomainVoltage('5 * x(t)').fourier(), FourierDomainVoltage('5 * X(f)'), "5 * x(t)")

    def test_inverse_fourier(self):

        self.assertEqual((f * 0 + 1).inverse_fourier(), DiracDelta(t),
                         "1")
        self.assertEqual((f * 0 + 10).inverse_fourier(), 10
                         * DiracDelta(t), "0")
        self.assertEqual(FourierDomainVoltage('V(f)').inverse_fourier(), TimeDomainVoltage('v(t)'), "V(f)")
        self.assertEqual(FourierDomainVoltage('10 * V(f)').inverse_fourier(),
                         TimeDomainVoltage('10 * v(t)'), "V(f)")
        self.assertEqual(TimeDomainVoltage('v(t)').fourier().inverse_fourier(),
                         TimeDomainVoltage('v(t)'), "v(t)")
        self.assertEqual(TimeDomainVoltage('v(t/2)').fourier().inverse_fourier(),
                         TimeDomainVoltage('v(t/2)'), "v(t/2)")        
                         
    def test_fourier2(self):

        self.assertEqual((t * 0 + 1).fourier(), DiracDelta(f))
        self.assertEqual((t * 0 + 1).fourier().inverse_fourier(), 1)
        self.assertEqual(t.fourier(), 2 * j * pi * f * DiracDelta(f, 1))
        self.assertEqual(2 * cos(2 * pi * t).fourier(),
                         DiracDelta(f - 1) + DiracDelta(f + 1))
        self.assertEqual(2 * sin(2 * pi * t).fourier(),
                         -j * DiracDelta(f - 1) + j * DiracDelta(f + 1))
        self.assertEqual(exp(j * 2 * pi * t).fourier(), DiracDelta(f - 1))
        self.assertEqual(exp(j * 2 * pi * t).fourier().inverse_fourier(),
                         exp(j * 2 * pi * t))
        self.assertEqual(exp(-j * 2 * pi * t).fourier(), DiracDelta(f + 1))
        self.assertEqual(exp(-j * 2 * pi * t).fourier().inverse_fourier(),
                         exp(-j * 2 * pi * t))

    def test_inverse_fourier2(self):

        self.assertEqual((1 / (s + 1))(j * omega, causal=True).inverse_fourier(), exp(-t) * Heaviside(t))
        self.assertEqual((1 / (s + 1))(j * omega, causal=True)(2 * pi * f).inverse_fourier(), exp(-t) * Heaviside(t))

