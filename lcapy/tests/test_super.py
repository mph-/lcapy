from lcapy import *
from lcapy.voltage import Vname
from lcapy.current import Iname
from lcapy.superpositionvoltage import SuperpositionVoltage
from lcapy.superpositioncurrent import SuperpositionCurrent
from lcapy.phasor import PhasorDomainExpression
import unittest
import sympy as sym


def div(a, b):
    return a / b


def mul(a, b):
    return a * b


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def assertEqual2(self, ans1, ans2, comment):

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            ans1.pprint()
            ans2.pprint()
            raise AssertionError(e)

    def test_voltage_properties(self):
        self.assertEqual(SuperpositionVoltage(
            3).is_dc, True, "Voltage(3).is_dc")
        self.assertEqual(SuperpositionVoltage(phasorvoltage(
            3)).is_ac, True, "Voltage(Vphasor(3)).is_ac")
        self.assertEqual(SuperpositionVoltage(voltage(2), phasorvoltage(3)).is_ac, False,
                         "Voltage(Vconst(2), Vphasor(3)).is_ac")
        self.assertEqual(SuperpositionVoltage(voltage(2), phasorvoltage(3)).is_ac, False,
                         "Voltage(Vconst(2), Vphasor(3)).is_dc")

    def test_voltage_add_sub_dc(self):
        self.assertEqual2(SuperpositionVoltage(
            3).dc, voltage(3), "Voltage(3).dc")
        self.assertEqual2(SuperpositionVoltage(2, 3).dc,
                          voltage(5), "Voltage(2, 3).dc")
        self.assertEqual2(SuperpositionVoltage(
            2, 3).ac, {}, "Voltage(2, 3).ac")
        self.assertEqual2(-SuperpositionVoltage(2).dc,
                          voltage(-2), "-Voltage(2).dc")
        self.assertEqual2(SuperpositionVoltage(2) + SuperpositionVoltage(3), SuperpositionVoltage(5),
                          "Voltage(2) + Voltage(3)")
        self.assertEqual2(SuperpositionVoltage(2) - SuperpositionVoltage(3), SuperpositionVoltage(-1),
                          "Voltage(2) - Voltage(3)")

    def test_current_add_sub_dc(self):
        self.assertEqual2(SuperpositionCurrent(
            3).dc, current(3), "Current(3).dc")
        self.assertEqual2(SuperpositionCurrent(2, 3).dc,
                          current(5), "Current(2, 3).dc")
        self.assertEqual2(SuperpositionCurrent(
            2, 3).ac, {}, "Current(2, 3).ac")
        self.assertEqual2(-SuperpositionCurrent(2).dc,
                          current(-2), "-Current(2).dc")
        self.assertEqual2(SuperpositionCurrent(2) + SuperpositionCurrent(3), SuperpositionCurrent(5),
                          "Current(2) + Current(3)")
        self.assertEqual2(SuperpositionCurrent(2) - SuperpositionCurrent(3), SuperpositionCurrent(-1),
                          "Current(2) - Current(3)")

    # def test_Voltage_noise(self):
        # self.assertEqual((AngularFourierNoiseDomainVoltage(3) + AngularFourierNoiseDomainVoltage(4)).expr, AngularFourierNoiseDomainVoltage(5).expr, "Vnoisy(3) + Vnoisy(4)")
        # self.assertEqual((SuperpositionVoltage(AngularFourierNoiseDomainVoltage(3)) + SuperpositionVoltage(AngularFourierNoiseDomainVoltage(4))).n.expr,
        #                   SuperpositionVoltage(AngularFourierNoiseDomainVoltage(5)).n.expr,
        #                   "Voltage(Vnoisy(3)) + Voltage(Vnoisy(4))")

    def test_voltage_has(self):

        a = SuperpositionVoltage('3 * exp(-t) * t * a')
        self.assertEqual(a.has(3), True, "has(3)")
        self.assertEqual(a.has(4), False, "has(4)")
        self.assertEqual(a.has(t), True, "has(t)")
        self.assertEqual(a.has_symbol(t), True, "has_symbol(t)")
        self.assertEqual(a.has_symbol('a'), True, "has_symbol(a)")
        self.assertEqual(a.has_symbol('b'), False, "has_symbol(b)")

    def test_voltage_transform(self):

        V1 = SuperpositionVoltage('3 * exp(-2 * t)')
        self.assertEqual(V1.transform(s), voltage(3 / (s + 2)), 'transform(s)')

        V2 = SuperpositionVoltage('3 * exp(-2 * t) * u(t)')
        self.assertEqual(V2.transform(s), voltage(3 / (s + 2)), 'transform(s)')
        self.assertEqual(simplify(V2.transform(
            f) - voltage(3 / (j * 2 * pi * f + 2))), 0, 'transform(f)')

    def test_voltage_subs(self):

        a = SuperpositionVoltage('V1')
        b = a.subs('V1', 1)
        c = SuperpositionVoltage(1)
        self.assertEqual(b, c, "Voltage.subs")

    def test_voltage_decompose(self):

        V1 = SuperpositionVoltage('1 + 3 * u(t) + cos(2 * pi * 3 * t)')
        self.assertEqual(V1.dc, voltage(1), '.dc')
        self.assertEqual(V1.transient, voltage('3 * u(t)'), '.transient')

    def test_voltage_oneport(self):

        V1 = V(3)

        self.assertEqual(V1.V.oneport().V, V1.V, 'oneport')

    def test_current_oneport(self):

        I1 = I(3)

        self.assertEqual(I1.I.oneport().I, I1.I, 'oneport')

    def test_Vname(self):

        self.assertEqual(Vname('V', 't'), voltage('v(t)'), 'v(t)')
        self.assertEqual(Vname('V', 's'), voltage('V(s)'), 'V(s)')
        # TODO: remove cache requirement
        self.assertEqual(Vname('V', 'dc', cache=True), voltage('V'), 'V')

    def test_Iname(self):

        self.assertEqual(Iname('I', 't'), current('i(t)'), 'i(t)')
        self.assertEqual(Iname('I', 's'), current('I(s)'), 'I(s)')
        self.assertEqual(Iname('I', 'dc', cache=True), current('I'), 'I')

    def test_voltage_phasor(self):

        V = SuperpositionVoltage(3 * sin(7 * t) + 2 * cos(14 * t))
        self.assertEqual(V[7].magnitude, voltage(3), 'magnitude')
        self.assertEqual(V[14].omega, 14, 'omega')

    def test_super_printing(self):

        self.assertEqual(SuperpositionVoltage(3).latex(),
                         '3', 'SuperpositionVoltage(3)')
        self.assertEqual(SuperpositionVoltage(3 * s).latex(),
                         '3 s', 'SuperpositionVoltage(3 * s)')
        self.assertEqual(SuperpositionVoltage(2 * sin(omega * t)).latex(),
                         '2 \\sin{\\left(\\omega t \\right)}',
                         'SuperpositionVoltage(2 * sin(omega * t))')
        self.assertEqual(SuperpositionVoltage(3 + 2 * sin(omega * t)).latex(),
                         '2 \\sin{\\left(\\omega t \\right)} + 3',
                         'SuperpositionVoltage(3 + 2 * sin(omega * t))')

    def test_super_op_const(self):

        V = SuperpositionVoltage(10)
        Z = impedance(2)
        Y = admittance(1 / 2)
        I = SuperpositionCurrent(5)

        self.assertEqual(V._div(Z), I, 'V / R')
        self.assertEqual(I._div(Y), V, 'I / G')
        self.assertEqual(V._mul(Y), I, 'V / G')
        self.assertEqual(I._mul(Z), V, 'I * R')

    def test_super_op_laplace(self):

        V = SuperpositionVoltage(10 / s)
        Z = impedance(2 * s)
        Y = admittance(1 / (2 * s))
        I = SuperpositionCurrent(5 / s**2)

        self.assertEqual(V._div(Z), I, 'V / R')
        self.assertEqual(I._div(Y), V, 'I / G')
        self.assertEqual(V._mul(Y), I, 'V / G')
        self.assertEqual(I._mul(Z), V, 'I * R')

    def test_super_op_phasor(self):

        V = SuperpositionVoltage((10 * sin(omega * t)).as_phasor())
        # TODO, try with jw
        Z = impedance(2 * s)
        Y = admittance(1 / (2 * s))
        I = SuperpositionCurrent((-5 / omega * cos(omega * t)).as_phasor())

        self.assertEqual(V._div(Z), I, 'V / R')
        self.assertEqual(I._div(Y), V, 'I / G')
        self.assertEqual(V._mul(Y), I, 'V / G')
        self.assertEqual(I._mul(Z), V, 'I * R')

    def test_super_time_div_laplace(self):

        self.assertRaises(ValueError, div, voltage('v(t)'), impedance(3 * s))

    def test_super_time_mul_laplace(self):

        self.assertRaises(ValueError, div, current('i(t)'), admittance(3 * s))

    def test_super_mul(self):

        self.assertEqual(SuperpositionCurrent(5) * impedance(2),
                         SuperpositionVoltage(10), 'I * R')
        self.assertEqual(SuperpositionCurrent(
            5) * impedance(0 * s + 2), SuperpositionVoltage(10), 'I * R')

    def test_super_reverse_mul(self):

        self.assertEqual(impedance(2) * SuperpositionCurrent(5),
                         SuperpositionVoltage(10), 'I * R')

    def test_super_sympy_float(self):

        V = SuperpositionVoltage(sym.Float(8))

        self.assertEqual(isinstance(V.dc.expr, sym.Rational),
                         True, 'Float -> Rational')

    def test_super_phasor(self):

        V = SuperpositionVoltage(PhasorDomainExpression(3 + 4j))
        self.assertEqual(V.magnitude, 5, 'Magnitude of single phasor')
        self.assertEqual(V.phase, atan2(4, 3), 'Phase of single phasor')
