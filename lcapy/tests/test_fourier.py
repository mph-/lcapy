from lcapy import *
from lcapy.fexpr import FourierDomainExpression
from lcapy.expr import TimeDomainExpression
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_fourier(self):

        self.assertEqual(DiracDelta(t).fourier(), 1, "DiracDelta(t)")
        self.assertEqual(TimeDomainExpression('x(t)').fourier(), FourierDomainExpression('X(f)'), "x(t)")
        self.assertEqual(TimeDomainExpression('5 * x(t)').fourier(), FourierDomainExpression('5 * X(f)'), "5 * x(t)")

    def test_inverse_fourier(self):

        self.assertEqual((f * 0 + 1).inverse_fourier(), DiracDelta(t),
                         "1")
        self.assertEqual((f * 0 + 10).inverse_fourier(), 10
                         * DiracDelta(t), "0")
        self.assertEqual(FourierDomainExpression('V(f)').inverse_fourier(), TimeDomainExpression('v(t)'), "V(f)")
        self.assertEqual(FourierDomainExpression('10 * V(f)').inverse_fourier(),
                         TimeDomainExpression('10 * v(t)'), "V(f)")
        self.assertEqual(TimeDomainExpression('v(t)').fourier().inverse_fourier(),
                         TimeDomainExpression('v(t)'), "v(t)")
        self.assertEqual(TimeDomainExpression('v(t/2)').fourier().inverse_fourier(),
                         TimeDomainExpression('v(t/2)'), "v(t/2)")        
                         
    def test_fourier2(self):

        # 1
        self.assertEqual((t * 0 + 1).fourier(), DiracDelta(f))
        self.assertEqual((t * 0 + 1).fourier().inverse_fourier(), 1)

        # t
        self.assertEqual(t.fourier(), j * DiracDelta(f, 1) / (2 * pi))
        self.assertEqual(t.fourier().inverse_fourier(), t)

        # t**2
        self.assertEqual((t**2).fourier(), -DiracDelta(f, 2) / (2 * pi)**2)
        self.assertEqual((t**2).fourier().inverse_fourier(), t**2)
        
        # 2 * cos(2 * pi * t)
        self.assertEqual(2 * cos(2 * pi * t).fourier(),
                         DiracDelta(f - 1) + DiracDelta(f + 1))

        # 2 * sin(2 * pi * t)        
        self.assertEqual(2 * sin(2 * pi * t).fourier(),
                         -j * DiracDelta(f - 1) + j * DiracDelta(f + 1))

        # exp(j * 2 * pi * t)
        self.assertEqual(exp(j * 2 * pi * t).fourier(), DiracDelta(f - 1))
        self.assertEqual(exp(j * 2 * pi * t).fourier().inverse_fourier(),
                         exp(j * 2 * pi * t))

        # exp(-j * 2 * pi * t)        
        self.assertEqual(exp(-j * 2 * pi * t).fourier(), DiracDelta(f + 1))
        self.assertEqual(exp(-j * 2 * pi * t).fourier().inverse_fourier(),
                         exp(-j * 2 * pi * t))

        # u(t)
        self.assertEqual(u(t).fourier(), DiracDelta(f) / 2 + 1 / (j * 2 * pi * f))
        # TODO teach SymPy that sign(t) + 1 = 2 * u(t)
        #self.assertEqual(u(t).fourier().inverse_fourier(), u(t))

        # t * u(t)
        self.assertEqual((t * u(t)).fourier(), -1 / (2 * pi * f)**2 + j * DiracDelta(f, 1) / (4 * pi))
        # TODO teach SymPy that sign(t) + 1 = 2 * u(t)        
        #self.assertEqual((t * u(t)).fourier().inverse_fourier(), t * u(t))        

        # sign(t)
        self.assertEqual(sign(t).fourier(), 1 / (j * pi * f))
        self.assertEqual(sign(t).fourier().inverse_fourier(), sign(t))

        # t * sign(t)
        self.assertEqual((t * sign(t)).fourier(), -1 / (2 * pi**2 * f**2))
        self.assertEqual((t * sign(t)).fourier().inverse_fourier(), t * sign(t))

        # 2 * cos(2 * pi * t + 3)
        self.assertEqual(2 * cos(2 * pi * t + 3).fourier(),
                         DiracDelta(f - 1) * exp(3 * j) +
                         DiracDelta(f + 1) * exp(-3 * j))
        
