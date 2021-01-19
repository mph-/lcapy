from lcapy import *
import unittest
from sympy import Rational


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_LSection(self):
        """Lcapy: check LSection

        """
        a = LSection(R(10), R(30))

        self.assertEqual(a.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(a.Igain12, -1, "Igain12 incorrect.")
        self.assertEqual(a.Z1oc, impedance(40), "Z1oc incorrect.")
        self.assertEqual(a.Z2oc, impedance(30), "Z2oc incorrect.")
        self.assertEqual(a.Z1sc, impedance(10), "Z1sc incorrect.")
        self.assertEqual(a.Z2sc, impedance(7.5), "Z2sc incorrect.")
        self.assertEqual(a.Vgain(1, 2), 0.75, "Vgain 1->2 incorrect.")
        self.assertEqual(a.Vgain(2, 1), 1, "Vgain 2->1 incorrect.")
        self.assertEqual(a.Igain(1, 2), -1, "Igain 1->2 incorrect.")
        self.assertEqual(a.Igain(2, 1), -0.75, "Igain 2->1 incorrect.")
        self.assertEqual(
            a.Vresponse(Vdc(1), 1, 2), voltage(0.75 / s), "Vresponse 1->2 incorrect.")
        self.assertEqual(
            a.Vresponse(Vdc(1), 2, 1), voltage(1 / s), "Vresponse 2->1 incorrect.")
        self.assertEqual(
            a.Iresponse(Idc(1), 1, 2), current(-1 / s), "Iresponse 1->2 incorrect.")
        self.assertEqual(
            a.Iresponse(Idc(1), 2, 1), current(-0.75 / s), "Iresponse 2->1 incorrect.")

        b = a.load(R(30))
        self.assertEqual(b.Z, impedance(25), "R loaded Z incorrect.")
        self.assertEqual(b.Voc.dc, 0, "R loaded V incorrect.")

        c = a.load(R(30) + Vdc(5))
        self.assertEqual(c.Z, impedance(25), "T loaded Z incorrect.")
        self.assertEqual(c.Voc(s), voltage(2.5 / s), "T loaded V incorrect.")

        d = LSection(R(10), R(30) + Vdc(6))
        self.assertEqual(d.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(d.Igain12, -1, "Igain12 incorrect.")
        self.assertEqual(d.Z1oc, impedance(40), "Z1oc incorrect.")
        self.assertEqual(d.Z2oc, impedance(30), "Z2oc incorrect.")
        self.assertEqual(d.V1oc, voltage(6 / s), "V1oc incorrect.")
        self.assertEqual(d.V2oc, voltage(6 / s), "V2oc incorrect.")
        self.assertEqual(d.I1sc, 0 / s, "I2sc incorrect.")
        self.assertEqual(d.I2sc, current(-0.2 / s), "I2sc incorrect.")

        e = d.chain(Series(R(10)))
        self.assertEqual(e.Z1oc, impedance(40), "Z1oc incorrect.")
        self.assertEqual(e.Z2oc, impedance(40), "Z2oc incorrect.")
        self.assertEqual(e.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(e.Igain12, -0.75, "Igain12 incorrect.")

    def test_Shunt_parallel(self):
        """Lcapy: check Shunts in parallel

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.chain(Shunt(R(30)))
        self.assertEqual(b.Z1oc, impedance(7.5), "Z1oc incorrect.")
        self.assertEqual(b.Z2oc, impedance(7.5), "Z2oc incorrect.")
        self.assertEqual(b.V1oc, voltage(3.75 / s), "V1oc incorrect.")
        self.assertEqual(b.V2oc, voltage(3.75 / s), "V2oc incorrect.")

    def test_Series_series(self):
        """Lcapy: check Series in series

        """

        a = Series(R(10) + Vdc(5))
        b = a.chain(Series(R(30)))
        self.assertEqual(b.Y1sc, admittance(0.025), "Y1sc incorrect.")
        self.assertEqual(b.Y2sc, admittance(0.025), "Y2sc incorrect.")
        # This needs some thinking...
        # self.assertEqual(b.I1sc, -0.125 / s, "I1sc incorrect.")
        # self.assertEqual(b.I2sc, 0.125 / s, "I2sc incorrect.")

        c = a.chain(a)

    def test_LSection_models(self):
        """Lcapy: check LSection models

        """
        a = LSection(R(10) + Vdc(5), R(20))

        self.assertEqual(a.Bparams.B11, 1, "incorrect B11.")
        self.assertEqual(a.Bparams.B12, -10, "incorrect B12.")
        # Note, -1 / 20 = -0.5 cannot be represented exactly as a float
        self.assertEqual(a.Bparams.B21, Rational(-1) / 20, "incorrect B21.")
        self.assertEqual(a.Bparams.B22, 1.5, "incorrect B22.")

        b = a.Zmodel

        self.assertEqual(b.Bparams.B11, 1, "incorrect B11.")
        self.assertEqual(b.Bparams.B12, -10, "incorrect B12.")
        self.assertEqual(b.Bparams.B21, Rational(-1) / 20, "incorrect B21.")
        self.assertEqual(b.Bparams.B22, 1.5, "incorrect B22.")

        c = a.Ymodel

        self.assertEqual(c.Bparams.B11, 1, "incorrect B11.")
        self.assertEqual(c.Bparams.B12, -10, "incorrect B12.")
        self.assertEqual(c.Bparams.B21, Rational(-1) / 20, "incorrect B21.")
        self.assertEqual(c.Bparams.B22, 1.5, "incorrect B22.")

        d = a.Hmodel

        self.assertEqual(d.Bparams.B11, 1, "incorrect B11.")
        self.assertEqual(d.Bparams.B12, -10, "incorrect B12.")
        self.assertEqual(d.Bparams.B21, Rational(-1) / 20, "incorrect B21.")
        self.assertEqual(d.Bparams.B22, 1.5, "incorrect B22.")

    def test_load(self):
        """Lcapy: check load

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.load(R(30))

        self.assertEqual(b.Z, impedance(7.5), "Shunt loaded R incorrect Z.")

    def test_open_circuit(self):
        """Lcapy: check open_circuit

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.open_circuit()

        self.assertEqual(b.Z, impedance(10), "incorrect Z.")
        self.assertEqual(b.Voc(s), voltage(5 / s), "incorrect V.")

    def test_short_circuit(self):
        """Lcapy: check short_circuit

        """

        a = Series(R(10) + Vdc(5))
        b = a.short_circuit()

        self.assertEqual(b.Z, impedance(10), "incorrect Z.")
        self.assertEqual(b.Voc(s), voltage(5 / s), "incorrect V.")

    def test_Shunt(self):
        """Lcapys: check Shunt

        """
        a = Shunt(R(10) + Vdc(5))
        self.assertEqual(a.Z1oc, impedance(10), "Z1oc incorrect.")
        self.assertEqual(a.Z2oc, impedance(10), "Z2oc incorrect.")
        self.assertEqual(a.V2oc, voltage(5 / s), "V2oc incorrect.")
        self.assertEqual(a.V1oc, voltage(5 / s), "V1oc incorrect.")
        # Cannot determine I1sc, I2sc, Ymn

        b = a.chain(Shunt(R(30)))
        self.assertEqual(b.Z1oc, impedance(7.5), "Z1oc incorrect.")
        self.assertEqual(b.Z2oc, impedance(7.5), "Z2oc incorrect.")
        self.assertEqual(b.V2oc, voltage(3.75 / s), "V2oc incorrect.")
        self.assertEqual(b.V1oc, voltage(3.75 / s), "V1oc incorrect.")

        c = a.chain(Series(R(30)))
        self.assertEqual(c.Z1oc, impedance(10), "Z1oc incorrect.")
        self.assertEqual(c.Z2oc, impedance(40), "Z2oc incorrect.")
        self.assertEqual(c.V1oc, voltage(5 / s), "V1oc incorrect.")
        print(c.V2oc)
        self.assertEqual(c.V2oc, voltage(5 / s), "V2oc incorrect.")
        self.assertEqual(c.I1sc, current(-0.5 / s), "I1sc incorrect.")
        self.assertEqual(c.I2sc, current(0 / s), "I2sc incorrect.")

        d = Shunt(R(10)).chain(Series(R(30)))

    def test_Series(self):
        """Lcapys: check Series

        """
        a = Series(R(10) + Vdc(5))
        self.assertEqual(a.Y1sc, admittance(0.1), "Y1sc incorrect.")
        self.assertEqual(a.Y2sc, admittance(0.1), "Y2sc incorrect.")
        self.assertEqual(a.I1sc, current(0.5 / s), "I1sc incorrect.")
        self.assertEqual(a.I2sc, current(-0.5 / s), "I2sc incorrect.")
        # Cannot determine V1oc, V2oc, Zmn

    def test_LSection(self):
        """Lcapys: check LSection

        """
        a = LSection(R(10), R(30))

        self.assertEqual(a.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(a.Igain12, -1, "Igain12 incorrect.")
        self.assertEqual(a.Z1oc, impedance(40), "Z1oc incorrect.")
        self.assertEqual(a.Z2oc, impedance(30), "Z2oc incorrect.")
        self.assertEqual(a.Z1sc, impedance(10), "Z1sc incorrect.")
        self.assertEqual(a.Z2sc, impedance(7.5), "Z2sc incorrect.")
        self.assertEqual(a.Vgain(1, 2), 0.75, "Vgain 1->2 incorrect.")
        self.assertEqual(a.Vgain(2, 1), 1, "Vgain 2->1 incorrect.")
        self.assertEqual(a.Igain(1, 2), -1, "Igain 1->2 incorrect.")
        self.assertEqual(a.Igain(2, 1), -0.75, "Igain 2->1 incorrect.")
        self.assertEqual(
            a.Vresponse(Vdc(1), 1, 2), voltage(0.75 / s), "Vresponse 1->2 incorrect.")
        self.assertEqual(
            a.Vresponse(Vdc(1), 2, 1), voltage(1 / s), "Vresponse 2->1 incorrect.")
        self.assertEqual(
            a.Iresponse(Idc(1), 1, 2), current(-1 / s), "Iresponse 1->2 incorrect.")
        self.assertEqual(
            a.Iresponse(Idc(1), 2, 1), current(-0.75 / s), "Iresponse 2->1 incorrect.")

        b = a.load(R(30))
        self.assertEqual(b.Z, impedance(25), "R loaded Z incorrect.")
        self.assertEqual(b.Voc.dc, 0, "R loaded V incorrect.")

        c = a.load(R(30) + Vdc(5))
        self.assertEqual(c.Z, impedance(25), "T loaded Z incorrect.")
        self.assertEqual(c.Voc(s), voltage(2.5 / s), "T loaded V incorrect.")

        d = LSection(R(10), R(30) + Vdc(6))
        self.assertEqual(d.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(d.Igain12, -1, "Igain12 incorrect.")
        self.assertEqual(d.Z1oc, impedance(40), "Z1oc incorrect.")
        self.assertEqual(d.Z2oc, impedance(30), "Z2oc incorrect.")
        self.assertEqual(d.V1oc, voltage(6 / s), "V1oc incorrect.")
        self.assertEqual(d.V2oc, voltage(6 / s), "V2oc incorrect.")
        self.assertEqual(d.I1sc, current(0 / s), "I2sc incorrect.")
        self.assertEqual(d.I2sc, current(-0.2 / s), "I2sc incorrect.")

        e = d.chain(Series(R(10)))
        self.assertEqual(e.Z1oc, impedance(40), "Z1oc incorrect.")
        self.assertEqual(e.Z2oc, impedance(40), "Z2oc incorrect.")
        self.assertEqual(e.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(e.Igain12, -0.75, "Igain12 incorrect.")

    def test_Shunt_parallel(self):
        """Lcapys: check Shunts in parallel

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.chain(Shunt(R(30)))
        self.assertEqual(b.Z1oc, impedance(7.5), "Z1oc incorrect.")
        self.assertEqual(b.Z2oc, impedance(7.5), "Z2oc incorrect.")
        self.assertEqual(b.V1oc, voltage(3.75 / s), "V1oc incorrect.")
        self.assertEqual(b.V2oc, voltage(3.75 / s), "V2oc incorrect.")

    def test_Series_series(self):
        """Lcapys: check Series in series

        """

        a = Series(R(10) + Vdc(5))
        b = a.chain(Series(R(30)))
        self.assertEqual(b.Y1sc, admittance(0.025), "Y1sc incorrect.")
        self.assertEqual(b.Y2sc, admittance(0.025), "Y2sc incorrect.")
        # The following fail and need some thinking...
        # self.assertEqual(b.I1sc, -0.125 / s, "I1sc incorrect.")
        # self.assertEqual(b.I2sc, 0.125 / s, "I2sc incorrect.")

    def test_load(self):
        """Lcapys: check load

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.load(R(30))

        self.assertEqual(b.Z, impedance(7.5), "Shunt loaded R incorrect Z.")

    def test_open_circuit(self):
        """Lcapys: check open_circuit

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.open_circuit()

        self.assertEqual(b.Z, impedance(10), "incorrect Z.")
        self.assertEqual(b.Voc(s), voltage(5 / s), "incorrect V.")

    def test_short_circuit(self):
        """Lcapys: check short_circuit

        """

        a = Series(R(10) + Vdc(5))
        b = a.short_circuit()

        self.assertEqual(b.Z, impedance(10), "incorrect Z.")
        self.assertEqual(b.Voc(s), voltage(5 / s), "incorrect V.")

    def test_LSection_models(self):
        """Lcapys: check LSection models

        """
        a = LSection(R(10) + Vdc(5), R(20))

        self.assertEqual(a.Bparams.B11, 1, "incorrect B11.")
        self.assertEqual(a.Bparams.B12, -10, "incorrect B12.")
        self.assertEqual(a.Bparams.B21, Rational(-1) / 20, "incorrect B21.")
        self.assertEqual(a.Bparams.B22, 1.5, "incorrect B22.")

        b = a.Zmodel

        self.assertEqual(b.Bparams.B11, 1, "incorrect B11.")
        self.assertEqual(b.Bparams.B12, -10, "incorrect B12.")
        self.assertEqual(b.Bparams.B21, Rational(-1) / 20, "incorrect B21.")
        self.assertEqual(b.Bparams.B22, 1.5, "incorrect B22.")

        c = a.Ymodel

        self.assertEqual(c.Bparams.B11, 1, "incorrect B11.")
        self.assertEqual(c.Bparams.B12, -10, "incorrect B12.")
        self.assertEqual(c.Bparams.B21, Rational(-1) / 20, "incorrect B21.")
        self.assertEqual(c.Bparams.B22, 1.5, "incorrect B22.")

        d = a.Hmodel

        self.assertEqual(d.Bparams.B11, 1, "incorrect B11.")
        self.assertEqual(d.Bparams.B12, -10, "incorrect B12.")
        self.assertEqual(d.Bparams.B21, Rational(-1) / 20, "incorrect B21.")
        self.assertEqual(d.Bparams.B22, 1.5, "incorrect B22.")

    def test_transforms(self):

        A = AMatrix.generic()        
        self.assertEqual(A.Aparams.Aparams.simplify(), A, "A.Aparams.Aparams")
        self.assertEqual(A.Bparams.Aparams.simplify(), A, "A.Bparams.Aparams")
        self.assertEqual(A.Gparams.Aparams.simplify(), A, "A.Gparams.Aparams")
        self.assertEqual(A.Hparams.Aparams.simplify(), A, "A.Hparams.Aparams")
        self.assertEqual(A.Sparams.Aparams.simplify(), A, "A.Sparams.Aparams")
        self.assertEqual(A.Tparams.Aparams.simplify(), A, "A.Tparams.Aparams")
        self.assertEqual(A.Yparams.Aparams.simplify(), A, "A.Yparams.Aparams")
        self.assertEqual(A.Zparams.Aparams.simplify(), A, "A.Zparams.Aparams")

        Z = ZMatrix.generic()            
        self.assertEqual(Z.Aparams.Zparams.simplify(), Z, "Z.Aparams.Zparams")
        self.assertEqual(Z.Bparams.Zparams.simplify(), Z, "Z.Bparams.Zparams")
        self.assertEqual(Z.Gparams.Zparams.simplify(), Z, "Z.Gparams.Zparams")
        self.assertEqual(Z.Hparams.Zparams.simplify(), Z, "Z.Hparams.Zparams")
        self.assertEqual(Z.Sparams.Zparams.simplify(), Z, "Z.Sparams.Zparams")
        self.assertEqual(Z.Tparams.Zparams.simplify(), Z, "Z.Tparams.Zparams")
        self.assertEqual(Z.Yparams.Zparams.simplify(), Z, "Z.Yparams.Zparams")
        self.assertEqual(Z.Zparams.Zparams.simplify(), Z, "Z.Zparams.Zparams")        
        
        
