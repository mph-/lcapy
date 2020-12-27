from lcapy import *
from lcapy.cexpr import ConstantExpression
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

    def test_tdomain_mul(self):

        v = voltage(5 * t)
        self.assertEqual(v * 3, 15 * t, "v(t) * const")
        self.assertEqual(3 * v, 15 * t, "const * v(t)")

    def test_sdomain_mul(self):
        
        V = voltage(5 / s)
        I = current(sexpr(5))        
        Y = admittance(s)
        Z = impedance(1 / s)        
        H = transfer(sexpr(10))
        U = sexpr(2)
        
        self.assertEqual(V * 3, 15 / s, "V(s) * const")
        self.assertEqual(3 * V, 15 / s, "const * V(s)")
        self.assertEqual(Y * V, I, "Y(s) * V(s)")
        self.assertEqual(V * Y, I, "V(s) * Y(s)")
        self.assertEqual(I * Z, V, "I(s) * Z(s)")
        self.assertEqual(Z * I, V, "Z(s) * I(s)")                
        self.assertEqual(V * H, voltage(50 / s), "V(s) * H(s)")
        self.assertEqual(H * V, voltage(50 / s), "H(s) * V(s)")
        self.assertEqual(I * H, current(sexpr(50)), "I(s) * H(s)")
        self.assertEqual(H * I, current(sexpr(50)), "H(s) * I(s)")
        self.assertEqual(U * U, sexpr(4), "U(s) * U(s)")        

    def test_fdomain_mul(self):
        
        V = voltage(5 / f)
        I = current(fexpr(5))        
        Y = admittance(f)
        Z = impedance(1 / f)                
        H = transfer(fexpr(10))
        U = fexpr(2)        
        
        self.assertEqual(V * 3, 15 / f, "V(f) * const")
        self.assertEqual(3 * V, 15 / f, "const * V(f)")
        self.assertEqual(Y * V, I, "Y(f) * V(f)")
        self.assertEqual(V * Y, I, "V(f) * Y(f)")
        self.assertEqual(I * Z, V, "I(f) * Z(f)")
        self.assertEqual(Z * I, V, "Z(f) * I(f)")                        
        self.assertEqual(V * H, voltage(50 / f), "V(f) * H(f)")
        self.assertEqual(H * V, voltage(50 / f), "H(f) * V(f)")
        self.assertEqual(I * H, current(fexpr(50)), "I(f) * H(f)")
        self.assertEqual(H * I, current(fexpr(50)), "H(f) * I(f)")
        self.assertEqual(U * U, fexpr(4), "U(f) * U(f)")                

    def test_sdomain_div(self):
        
        V = voltage(5 / s)
        I = current(sexpr(5))        
        Y = admittance(s)
        Z = impedance(1 / s)        
        H = transfer(sexpr(10))
        U = sexpr(2)        
        
        self.assertEqual(V / 5, 1 / s, "V(s) / 5")
        self.assertEqual(V / Z, I, "V(s) / Z(s)")
        self.assertEqual(I / Y, V, "I(s) / Y(s)")
        self.assertEqual(I / H, current(sexpr(1 / 2)), "I(s) / H(s)")
        self.assertEqual(V / H, voltage(sexpr(1 / 2 / s)), "V(s) / H(s)")
        self.assertEqual(U / U, sexpr(1), "U(s) / U(s)")        

    def test_sdomain_rdiv(self):
        
        Y = admittance(s)
        Z = impedance(1 / s)        
        H = transfer(sexpr(10))
        
        self.assertEqual(1 / Z, Y, "1 / Z(s)")
        self.assertEqual(1 / Y, Z, "1 / Y(s)") 
        
    def test_rmul(self):

        self.assertEqual2((j * omega).__class__,  AngularFourierDomainExpression, "j * omega fail.")

    def test_comparison(self):

        a = expr('3')
        b = expr('4')
        self.assertEqual(a < b, True, "a < b")
        self.assertEqual(a > b, False, "a > b")
        self.assertEqual(a <= b, True, "a <= b")
        self.assertEqual(a >= b, False, "a >= b")        
        self.assertEqual(a == b, False, "a == b")
        self.assertEqual(a != b, True, "a != b")        

    def test_rdiv(self):

        self.assertEqual(3 / expr('3'), 1, "3 / expr('3')")

    def test_div(self):

        self.assertEqual(expr('3') / 3, 1, "expr('3') / 3")        
        
    def test_rsub(self):

        self.assertEqual(3 - expr('3'), 0, "3 - expr('3')")

    def test_sub(self):

        self.assertEqual(expr('3') - 3, 0, "expr('3') - 3")

        
