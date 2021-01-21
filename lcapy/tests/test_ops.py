from lcapy import *
from lcapy.omegaexpr import AngularFourierDomainExpression
from lcapy.texpr import TimeDomainExpression
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
        self.assertEqual(v * 3, voltage(15 * t), "v(t) * const")
        self.assertEqual(3 * v, voltage(15 * t), "const * v(t)")

    def test_sdomain_mul(self):
        
        V = voltage(5 / s)
        I = current(sexpr(5))        
        Y = admittance(s)
        Z = impedance(1 / s)        
        H = transfer(sexpr(10))
        U = sexpr(2)
        C = cexpr(3)
        
        self.assertEqual(V * 3, voltage(15 / s), "V(s) * const")
        self.assertEqual(V * C, voltage(15 / s), "V(s) * const")        
        self.assertEqual(3 * V, voltage(15 / s), "const * V(s)")
        self.assertEqual(C * V, voltage(15 / s), "const * V(s)")        
        self.assertEqual(Y * V, I, "Y(s) * V(s)")
        self.assertEqual(V * Y, I, "V(s) * Y(s)")
        self.assertEqual(I * Z, V, "I(s) * Z(s)")
        self.assertEqual(Z * I, V, "Z(s) * I(s)")
        self.assertEqual(Z * Y, transfer(1), "Z(s) * Y(s)")
        self.assertEqual(Y * Z, transfer(1), "Y(s) * Z(s)")
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
        C = cexpr(3)
        
        self.assertEqual(V * 3, voltage(15 / f), "V(f) * const")
        self.assertEqual(V * C, voltage(15 / f), "V(f) * const")
        self.assertEqual(3 * V, voltage(15 / f), "const * V(f)")
        self.assertEqual(C * V, voltage(15 / f), "const * V(f)")        
        self.assertEqual(Y * V, I, "Y(f) * V(f)")
        self.assertEqual(V * Y, I, "V(f) * Y(f)")
        self.assertEqual(I * Z, V, "I(f) * Z(f)")
        self.assertEqual(Z * I, V, "Z(f) * I(f)")                        
        self.assertEqual(V * H, voltage(50 / f), "V(f) * H(f)")
        self.assertEqual(H * V, voltage(50 / f), "H(f) * V(f)")
        self.assertEqual(I * H, current(fexpr(50)), "I(f) * H(f)")
        self.assertEqual(H * I, current(fexpr(50)), "H(f) * I(f)")
        self.assertEqual(U * U, fexpr(4), "U(f) * U(f)")

    def test_mixed_mul(self):
        
        self.assertEqual(type(omega * t), TimeDomainExpression, "omega * t")
        self.assertEqual(type(t * omega), TimeDomainExpression, "omega * t")

    def test_sdomain_div(self):
        
        V = voltage(5 / s)
        I = current(sexpr(5))        
        Y = admittance(s)
        Z = impedance(1 / s)        
        H = transfer(sexpr(10))
        U = sexpr(2)        
        
        self.assertEqual(V / 5, voltage(1 / s), "V(s) / 5")
        self.assertEqual(V / Z, I, "V(s) / Z(s)")
        self.assertEqual(I / Y, V, "I(s) / Y(s)")
        self.assertEqual(I / H, current(sexpr(1 / 2)), "I(s) / H(s)")
        self.assertEqual(V / H, voltage(sexpr(1 / 2 / s)), "V(s) / H(s)")
        self.assertEqual(U / U, sexpr(1), "U(s) / U(s)")
        self.assertEqual(Z / Z, sexpr(1), "Z(s) / Z(s)")
        self.assertEqual(Y / Y, sexpr(1), "Y(s) / Y(s)")                        

    def test_sdomain_rdiv(self):
        
        Y = admittance(s)
        Z = impedance(1 / s)        
        H = transfer(sexpr(10))
        C = cexpr(1)        
        
        self.assertEqual(1 / Z, Y, "1 / Z(s)")
        self.assertEqual(1 / Y, Z, "1 / Y(s)")
        self.assertEqual(C / Z, Y, "C / Z(s)")
        self.assertEqual(C / Y, Z, "C / Y(s)")         

    def test_sdomain_add(self):
        
        V1 = voltage(s)
        V2 = voltage(s + 4)
        C = cexpr(3)
        Z = impedance(3 * s)
        CZ = impedance(3)        
        
        self.assertEqual(V1 + V2, voltage(2 * s + 4), "V1(s) + V2(s)")
        self.assertEqual(CZ + Z, impedance(3 * s + 3), "CZ + Z(s)")
        self.assertEqual(Z + CZ, impedance(3 * s + 3), "Z(s) + CZ")        

        loose = state.loose_units
        state.loose_units = True
        self.assertEqual(C + Z, impedance(3 * s + 3), "C + Z(s)")
        self.assertEqual(Z + C, impedance(3 * s + 3), "Z(s) + C")
        self.assertEqual(C + V1, voltage(3 + s), "C + V1(s)")
        self.assertEqual(V1 + C, voltage(3 + s), "V1(s) + C")
        state.loose_units = loose
        
    def test_fdomain_add(self):
        
        V1 = voltage(f)
        V2 = voltage(f + 4)
        C = cexpr(3)
        
        self.assertEqual(V1 + V2, voltage(2 * f + 4), "V1(f) + V2(f)")

        loose = state.loose_units        
        state.loose_units = True        
        self.assertEqual(C + V1, voltage(3 + f), "C + V1(f)")
        self.assertEqual(V1 + C, voltage(3 + f), "V1(f) + C")
        self.assertEqual(V1 + 3, voltage(3 + f), "V1(f) + 3")
        self.assertEqual(3 + V1, voltage(3 + f), "3 + V1(f)")
        state.loose_units = loose
        
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

        a = voltage(3)
        b = current(3)
        self.assertEqual(a != b, True, "a != b")

        a = voltage(3 * s)
        b = voltage(3)
        self.assertEqual(a != b, True, "a != b")

        self.assertEqual(jw != s, True, "jw != s")                

    def test_rdiv(self):

        self.assertEqual(3 / expr('3'), 1, "3 / expr('3')")

    def test_div(self):

        self.assertEqual(expr('3') / 3, 1, "expr('3') / 3")        
        
    def test_rsub(self):

        self.assertEqual(3 - expr('3'), 0, "3 - expr('3')")

    def test_sub(self):

        self.assertEqual(expr('3') - 3, 0, "expr('3') - 3")

    def test_time_op_const(self):

        self.assertEqual(voltage(texpr(10)) / impedance(2), current(5), 'V / R')
        self.assertEqual(current(texpr(5)) / admittance(1 / 2), voltage(10), 'V / G')
        self.assertEqual(voltage(texpr(10)) * admittance(1 / 2), current(5), 'V * Y')
        self.assertEqual(current(texpr(5)) * impedance(2), voltage(10), 'I * Z')

        # Invoke conversion of LaplaceDomain to ConstantDomain
        self.assertEqual(current(texpr(5)) * impedance(0 * s + 2), voltage(10), 'I * Z')
        self.assertEqual(voltage(texpr(10)) / impedance(0 * s + 2), current(5), 'V / Z')                
        

        
