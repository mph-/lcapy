from msignal.mcircuit import *
import unittest
import sympy as sym
from msignal.mrf import MRF

s = sym.var('s')

class McircuitsTester(unittest.TestCase):
    """Unit tests for mcircuits

    """

    def test_R(self):
        """Mcircuits: check R

        """
        a = R(10)
        self.assertEqual(a.Z, 10, "Z incorrect.")


    def test_L(self):
        """Mcircuits: check L

        """
        a = L(10, 5)
        self.assertEqual(a.Z, 10 * s, "Z incorrect.")
        self.assertEqual(a.V, -10 * 5, "V incorrect.")


    def test_C(self):
        """Mcircuits: check C

        """
        a = C(10, 5)

        self.assertEqual(a.Z, 1 / (10 * s), "Z incorrect.")
        self.assertEqual(a.V, 5 / s, "V incorrect.")


    def test_RplusL(self):
        """Mcircuits: check R + L

        """
        a = R(10) + L(5)
        self.assertEqual(a.Z, 5 * s + 10, "Z incorrect.")
        self.assertEqual(type(a), Thevenin, "type incorrect.")


    def test_Shunt(self):
        """Mcircuits: check Shunt

        """
        a = Shunt(R(10) + V(5))
        self.assertEqual(a.Z1oc, 10, "Z1oc incorrect.")
        self.assertEqual(a.Z2oc, 10, "Z2oc incorrect.")
        self.assertEqual(a.V2oc.mrf(), MRF(5, (1, 0)), "V2oc incorrect.")
        print(a.I2sc)
        self.assertEqual(a.I2sc.mrf(), MRF(0.5, (1, 0)), "I2sc incorrect.")
        self.assertEqual(a.V1oc.mrf(), MRF(5, (1, 0)), "V1oc incorrect.")
        self.assertEqual(a.I1sc.mrf(), MRF(0.5, (1, 0)), "I1sc incorrect.")

        b = a.chain(Shunt(R(30)))
        self.assertEqual(b.Z1oc.mrf(), MRF(7.5, 1), "Z1oc incorrect.")
        self.assertEqual(b.Z2oc.mrf(), MRF(7.5, 1), "Z2oc incorrect.")
        self.assertEqual(b.V2oc.mrf(), MRF(3.75, (1, 0)), "V2oc incorrect.")
        self.assertEqual(b.I2sc.mrf(), MRF(0.5, (1, 0)), "I2sc incorrect.")
        self.assertEqual(b.V1oc.mrf(), MRF(3.75, (1, 0)), "V1oc incorrect.")
        self.assertEqual(b.I1sc.mrf(), MRF(0.5, (1, 0)), "I1sc incorrect.")

        c = a.chain(Series(R(30)))
        self.assertEqual(c.Z1oc, 10, "Z1oc incorrect.")
        self.assertEqual(c.Z2oc.mrf(), MRF(40, 1), "Z2oc incorrect.")
        self.assertEqual(c.V1oc.mrf(), MRF(5, (1, 0)), "V1oc incorrect.")
        self.assertEqual(c.V2oc.mrf(), MRF(5, (1, 0)), "V2oc incorrect.")
        self.assertEqual(c.I1sc.mrf(), MRF(-0.5, (1, 0)), "I1sc incorrect.")
        self.assertEqual(c.I2sc.mrf(), MRF(0, (1, 0)), "I2sc incorrect.")

        d = Shunt(R(10)).chain(Series(R(30)))


    def test_Series(self):
        """Mcircuits: check Series

        """
        a = Series(R(10) + V(5))
        self.assertEqual(a.Y1sc, 0.1, "Y1sc incorrect.")
        self.assertEqual(a.Y1sc.mrf(), MRF(1, 10), "Y1sc incorrect.")
        self.assertEqual(a.Y2sc.mrf(), MRF(1, 10), "Y2sc incorrect.")
        print(a.V1oc)
        self.assertEqual(a.V1oc.mrf(), MRF(5, (1, 0)), "V1oc incorrect.")
        self.assertEqual(a.V2oc.mrf(), MRF(5, (1, 0)), "V2oc incorrect.")
        self.assertEqual(a.I1sc.mrf(), MRF(0.5, (1, 0)), "I1sc incorrect.")
        self.assertEqual(a.I2sc.mrf(), MRF(-0.5, (1, 0)), "I2sc incorrect.")


    def test_LSection(self):
        """Mcircuits: check LSection

        """
        a = LSection(R(10), R(30))

        self.assertEqual(a.Vgain12.mrf(), MRF(0.75, 1), "Vgain12 incorrect.")
        self.assertEqual(a.Igain12.mrf(), MRF(-1, 1), "Igain12 incorrect.")
        self.assertEqual(a.Z1oc.mrf(), MRF(40, 1), "Z1oc incorrect.")
        self.assertEqual(a.Z2oc.mrf(), MRF(30, 1), "Z2oc incorrect.")
        self.assertEqual(a.Z1sc, 10, "Z1sc incorrect.")
        self.assertEqual(a.Z2sc.mrf(), MRF(7.5, 1), "Z2sc incorrect.")
        self.assertEqual(a.Vgain(1, 2).mrf(), MRF(0.75, 1), "Vgain 1->2 incorrect.")
        self.assertEqual(a.Vgain(2, 1).mrf(), MRF(1, 1), "Vgain 2->1 incorrect.")
        self.assertEqual(a.Igain(1, 2).mrf(), MRF(-1, 1), "Igain 1->2 incorrect.")
        self.assertEqual(a.Igain(2, 1).mrf(), MRF(-0.75, 1), "Igain 2->1 incorrect.")
        self.assertEqual(a.Vresponse(V(1), 1, 2).mrf(), MRF(0.75, (1, 0)), "Vresponse 1->2 incorrect.")
        self.assertEqual(a.Vresponse(V(1), 2, 1).mrf(), MRF(1, (1, 0)), "Vresponse 2->1 incorrect.")
        self.assertEqual(a.Iresponse(I(1), 1, 2).mrf(), MRF(-1, (1, 0)), "Iresponse 1->2 incorrect.")
        self.assertEqual(a.Iresponse(I(1), 2, 1).mrf(), MRF(-0.75, (1, 0)), "Iresponse 2->1 incorrect.")

        b = a.load(R(30))
        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z.mrf(), MRF(25, 1), "R loaded Z incorrect.")
        self.assertEqual(b.V.mrf(), MRF(0, 1), "R loaded V incorrect.")

        c = a.load(R(30) + V(5))
        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(c.Z.mrf(), MRF(25, 1), "T loaded Z incorrect.")
        self.assertEqual(c.V.mrf(), MRF(2.5, (1, 0)), "T loaded V incorrect.")


        d = LSection(R(10), R(30) + V(6))
        self.assertEqual(d.Vgain12.mrf(), MRF(0.75, 1), "Vgain12 incorrect.")
        self.assertEqual(d.Igain12.mrf(), MRF(-1, 1), "Igain12 incorrect.")
        self.assertEqual(d.Z1oc.mrf(), MRF(40, 1), "Z1oc incorrect.")
        self.assertEqual(d.Z2oc.mrf(), MRF(30, 1), "Z2oc incorrect.")
        self.assertEqual(d.V1oc.mrf(), MRF(6, (1, 0)), "V1oc incorrect.")
        self.assertEqual(d.V2oc.mrf(), MRF(6, (1, 0)), "V2oc incorrect.")
        self.assertEqual(d.I1sc.mrf(), MRF(0, (1, 0)), "I2sc incorrect.")
        self.assertEqual(d.I2sc.mrf(), MRF(-0.2, (1, 0)), "I2sc incorrect.")
        

        e = d.chain(Series(R(10)))
        self.assertEqual(e.Z1oc.mrf(), MRF(40, 1), "Z1oc incorrect.")
        self.assertEqual(e.Z2oc.mrf(), MRF(40, 1), "Z2oc incorrect.")
        self.assertEqual(e.Vgain12.mrf(), MRF(0.75, 1), "Vgain12 incorrect.")
        self.assertEqual(e.Igain12.mrf(), MRF(-0.75, 1), "Igain12 incorrect.")


    def test_Shunt_parallel(self):
        """Mcircuits: check Shunts in parallel

        """

        a = Shunt(R(10) + V(5))
        b = a.chain(Shunt(R(30)))
        self.assertEqual(b.Z1oc.mrf(), MRF(7.5, 1), "Z1oc incorrect.")
        self.assertEqual(b.Z2oc.mrf(), MRF(7.5, 1), "Z2oc incorrect.")
        self.assertEqual(b.V1oc.mrf(), MRF(3.75, (1, 0)), "V1oc incorrect.")
        self.assertEqual(b.V2oc.mrf(), MRF(3.75, (1, 0)), "V2oc incorrect.")


    def test_Series_series(self):
        """Mcircuits: check Series in series

        """

        a = Series(R(10) + V(5))
        b = a.chain(Series(R(30)))
        self.assertEqual(b.Y1sc.mrf(), MRF(0.025, 1), "Y1sc incorrect.")
        self.assertEqual(b.Y2sc.mrf(), MRF(0.025, 1), "Y2sc incorrect.")
        self.assertEqual(b.I1sc.mrf(), MRF(-0.125, (1, 0)), "I1sc incorrect.")
        self.assertEqual(b.I2sc.mrf(), MRF(0.125, (1, 0)), "I2sc incorrect.")

        c = a.chain(a)



    def test_load(self):
        """Mcircuits: check load

        """

        a = Shunt(R(10) + V(5))
        b = a.load(R(30))

        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z.mrf(), MRF(7.5, 1), "Shunt loaded R incorrect Z.")
        

    def test_opencircuit(self):
        """Mcircuits: check opencircuit

        """

        a = Shunt(R(10) + V(5))
        b = a.opencircuit()

        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z, 10, "incorrect Z.")
        self.assertEqual(b.V.mrf(), MRF(5, (1, 0)), "incorrect V.")


    def test_shortcircuit(self):
        """Mcircuits: check shortcircuit

        """

        a = Series(R(10) + V(5))
        b = a.shortcircuit()

        self.assertEqual(type(b), Norton, "type incorrect.")
        self.assertEqual(b.Z, 10, "incorrect Z.")
        self.assertEqual(b.V.mrf(), MRF(5, (1, 0)), "incorrect V.")


    def test_LSection_models(self):
        """Mcircuits: check LSection models

        """
        a = LSection(R(10) + V(5), R(20))

        self.assertEqual(a.B.B11, 1, "incorrect B11.")
        self.assertEqual(a.B.B12, -10, "incorrect B12.")
        self.assertEqual(a.B.B21, -0.05, "incorrect B21.")
        self.assertEqual(a.B.B22, 1.5, "incorrect B22.")

        b = a.Zmodel

        self.assertEqual(b.B.B11, 1, "incorrect B11.")
        self.assertEqual(b.B.B12, -10, "incorrect B12.")
        self.assertEqual(b.B.B21, -0.05, "incorrect B21.")
        self.assertEqual(b.B.B22, 1.5, "incorrect B22.")

        c = a.Ymodel

        self.assertEqual(c.B.B11, 1, "incorrect B11.")
        self.assertEqual(c.B.B12, -10, "incorrect B12.")
        self.assertEqual(c.B.B21, -0.05, "incorrect B21.")
        self.assertEqual(c.B.B22, 1.5, "incorrect B22.")

        d = a.Hmodel

        self.assertEqual(d.B.B11, 1, "incorrect B11.")
        self.assertEqual(d.B.B12, -10, "incorrect B12.")
        self.assertEqual(d.B.B21, -0.05, "incorrect B21.")
        self.assertEqual(d.B.B22, 1.5, "incorrect B22.")
        
