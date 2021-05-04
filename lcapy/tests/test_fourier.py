from lcapy import *
from lcapy.fexpr import FourierDomainExpression
from lcapy.expr import TimeDomainExpression
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_fourier(self):

        self.assertEqual(DiracDelta(t).FT(), 1, "DiracDelta(t)")
        self.assertEqual(TimeDomainExpression('x(t)').FT(), FourierDomainExpression('X(f)'), "x(t)")
        self.assertEqual(TimeDomainExpression('5 * x(t)').FT(), FourierDomainExpression('5 * X(f)'), "5 * x(t)")

    def test_inverse_fourier(self):

        self.assertEqual((f * 0 + 1).IFT(), DiracDelta(t),
                         "1")
        self.assertEqual((f * 0 + 10).IFT(), 10
                         * DiracDelta(t), "0")
        self.assertEqual(FourierDomainExpression('V(f)').IFT(), TimeDomainExpression('v(t)'), "V(f)")
        self.assertEqual(FourierDomainExpression('10 * V(f)').IFT(),
                         TimeDomainExpression('10 * v(t)'), "V(f)")
        self.assertEqual(TimeDomainExpression('v(t)').FT().IFT(),
                         TimeDomainExpression('v(t)'), "v(t)")
        self.assertEqual(TimeDomainExpression('v(t/2)').FT().IFT(),
                         TimeDomainExpression('v(t/2)'), "v(t/2)")        
                         
    def test_fourier2(self):

        # 1
        self.assertEqual((t * 0 + 1).FT(), DiracDelta(f))
        self.assertEqual((t * 0 + 1).FT().IFT(), 1)

        # t
        self.assertEqual(t.FT(), j * DiracDelta(f, 1) / (2 * pi))
        self.assertEqual(t.FT().IFT(), t)

        # t**2
        self.assertEqual((t**2).FT(), -DiracDelta(f, 2) / (2 * pi)**2)
        self.assertEqual((t**2).FT().IFT(), t**2)
        
        # 2 * cos(2 * pi * t)
        self.assertEqual(2 * cos(2 * pi * t).FT(),
                         DiracDelta(f - 1) + DiracDelta(f + 1))

        # 2 * sin(2 * pi * t)        
        self.assertEqual(2 * sin(2 * pi * t).FT(),
                         -j * DiracDelta(f - 1) + j * DiracDelta(f + 1))

        # exp(j * 2 * pi * t)
        self.assertEqual(exp(j * 2 * pi * t).FT(), DiracDelta(f - 1))
        self.assertEqual(exp(j * 2 * pi * t).FT().IFT(),
                         exp(j * 2 * pi * t))

        # exp(-j * 2 * pi * t)        
        self.assertEqual(exp(-j * 2 * pi * t).FT(), DiracDelta(f + 1))
        self.assertEqual(exp(-j * 2 * pi * t).FT().IFT(),
                         exp(-j * 2 * pi * t))

        # u(t)
        self.assertEqual(u(t).FT(), DiracDelta(f) / 2 + 1 / (j * 2 * pi * f))
        # TODO teach SymPy that sign(t) + 1 = 2 * u(t)
        #self.assertEqual(u(t).FT().IFT(), u(t))

        # t * u(t)
        self.assertEqual((t * u(t)).FT(), -1 / (2 * pi * f)**2 + j * DiracDelta(f, 1) / (4 * pi))
        # TODO teach SymPy that sign(t) + 1 = 2 * u(t)        
        #self.assertEqual((t * u(t)).FT().IFT(), t * u(t))        

        # sign(t)
        self.assertEqual(sign(t).FT(), 1 / (j * pi * f))
        self.assertEqual(sign(t).FT().IFT(), sign(t))

        # t * sign(t)
        self.assertEqual((t * sign(t)).FT(), -1 / (2 * pi**2 * f**2))
        self.assertEqual((t * sign(t)).FT().IFT(), t * sign(t))

        # 2 * cos(2 * pi * t + 3)
        self.assertEqual(2 * cos(2 * pi * t + 3).FT(),
                         DiracDelta(f - 1) * exp(3 * j) +
                         DiracDelta(f + 1) * exp(-3 * j))

        self.assertEqual(2 * rect(t).FT(), 2 * sinc(f))
        self.assertEqual(rect(t / 2).FT(), 2 * sinc(2 * f))
        
        self.assertEqual(2 * sinc(t).FT(), 2 * rect(f))
        self.assertEqual(sinc(2 * t).FT(), rect(f / 2) / 2)

        self.assertEqual((rect(t) * cos(2 * pi * t)).FT(), sinc(f - 1) / 2 + sinc(f + 1) / 2)

        self.assertEqual(tri(t).FT(), sinc(f)**2)

        self.assertEqual(trap(t, 1).FT(), sinc(f)**2)
        self.assertEqual(trap(t, 0).FT(), sinc(f))
        self.assertEqual(trap(t, 0.5).FT(), 0.5 * sinc(0.5 * f) * sinc(f))
        
    def test_fourier_convolution(self):

        a = expr('Integral(3 * x(t - tau) * y(tau), (tau, -oo, oo))')
        A = a(f)
        self.assertEqual(A, expr('3 * X(f) * Y(f)'), "3 * X * Y") 
        #self.assertEqual(A(t), a, "a(f)(t)")       

        a = expr('Integral(3 * x(tau) * y(t - tau), (tau, -oo, oo))')
        A = a(f)        
        self.assertEqual(A, expr('3 * X(f) * Y(f)'), "3 * X * Y")
        #self.assertEqual(A(t), a, "a(f)(t)")               
