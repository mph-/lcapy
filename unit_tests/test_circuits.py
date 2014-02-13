from mcircuit import *
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


    def test_R_series_R(self):
        """Mcircuits: check R + R

        """
        a = R(10) + R(5)
        b = a.simplify()
        self.assertEqual(b.R, 15, "R incorrect.")
        self.assertEqual(b.Z, 15, "Z incorrect.")
        self.assertEqual(type(b), R, "type incorrect.")


    def test_L_series_L(self):
        """Mcircuits: check L + L

        """
        a = L(10) + L(5)
        b = a.simplify()
        self.assertEqual(b.L, 15, "L incorrect.")
        self.assertEqual(b.Z, 15 * s, "Z incorrect.")
        self.assertEqual(type(b), L, "type incorrect.")


    def test_C_series_C(self):
        """Mcircuits: check C + C

        """
        a = C(10) + C(15)
        b = a.simplify()
        self.assertEqual(b.C, 6, "C incorrect.")
        self.assertEqual(b.Z, 1 / (6 * s), "Z incorrect.")
        self.assertEqual(type(b), C, "type incorrect.")


    def test_V_series_V(self):
        """Mcircuits: check V + V

        """
        a = V(10) + V(5)
        b = a.simplify()
        self.assertEqual(b.v, 15, "V incorrect.")
        self.assertEqual(b.V, 15 / s, "V incorrect.")
        self.assertEqual(type(b), V, "type incorrect.")


    def test_R_series_L(self):
        """Mcircuits: check R + L

        """
        a = R(10) + L(5)
        self.assertEqual(a.Z, 5 * s + 10, "Z incorrect.")
        self.assertEqual(type(a), Thevenin, "type incorrect.")


    def test_R_parallel_R(self):
        """Mcircuits: check R | R

        """
        a = R(10) | R(15)
        b = a.simplify()
        self.assertEqual(b.R, 6, "R incorrect.")
        self.assertEqual(b.Z, 6, "Z incorrect.")
        self.assertEqual(type(b), R, "type incorrect.")


    def test_L_parallel_L(self):
        """Mcircuits: check L | L

        """
        a = L(10) | L(15)
        b = a.simplify()
        self.assertEqual(b.L, 6, "L incorrect.")
        self.assertEqual(b.Z, 6 * s, "Z incorrect.")
        self.assertEqual(type(b), L, "type incorrect.")


    def test_C_parallel_C(self):
        """Mcircuits: check C | C

        """
        a = C(10) | C(15)
        b = a.simplify()
        self.assertEqual(b.C, 25, "C incorrect.")
        self.assertEqual(b.Z, 1 / (25 * s), "Z incorrect.")
        self.assertEqual(type(b), C, "type incorrect.")


    def test_I_parallel_I(self):
        """Mcircuits: check I | I

        """
        a = I(10) | I(5)
        b = a.simplify()
        self.assertEqual(b.i, 15, "I incorrect.")
        self.assertEqual(b.I, 15 / s, "I incorrect.")
        self.assertEqual(type(b), I, "type incorrect.")


    def test_Shunt(self):
        """Mcircuits: check Shunt

        """
        a = Shunt(R(10) + V(5))
        self.assertEqual(a.Z1oc, 10, "Z1oc incorrect.")
        self.assertEqual(a.Z2oc, 10, "Z2oc incorrect.")
        self.assertEqual(a.V2oc, 5 / s, "V2oc incorrect.")
        print(a.I2sc)
        self.assertEqual(a.I2sc, 0.5 / s, "I2sc incorrect.")
        self.assertEqual(a.V1oc, 5 / s, "V1oc incorrect.")
        self.assertEqual(a.I1sc, 0.5 / s, "I1sc incorrect.")

        b = a.chain(Shunt(R(30)))
        self.assertEqual(b.Z1oc, 7.5, "Z1oc incorrect.")
        self.assertEqual(b.Z2oc, 7.5, "Z2oc incorrect.")
        self.assertEqual(b.V2oc, 3.75 / s, "V2oc incorrect.")
        self.assertEqual(b.I2sc, 0.5 / s, "I2sc incorrect.")
        self.assertEqual(b.V1oc, 3.75 / s, "V1oc incorrect.")
        self.assertEqual(b.I1sc, 0.5 / s, "I1sc incorrect.")

        c = a.chain(Series(R(30)))
        self.assertEqual(c.Z1oc, 10, "Z1oc incorrect.")
        self.assertEqual(c.Z2oc, 40, "Z2oc incorrect.")
        self.assertEqual(c.V1oc, 5 / s, "V1oc incorrect.")
        self.assertEqual(c.V2oc, 5 / s, "V2oc incorrect.")
        self.assertEqual(c.I1sc, -0.5 / s, "I1sc incorrect.")
        self.assertEqual(c.I2sc, 0 / s, "I2sc incorrect.")

        d = Shunt(R(10)).chain(Series(R(30)))


    def test_Series(self):
        """Mcircuits: check Series

        """
        a = Series(R(10) + V(5))
        self.assertEqual(a.Y1sc, 0.1, "Y1sc incorrect.")
        self.assertEqual(a.Y2sc, 0.1, "Y2sc incorrect.")
        print(a.V1oc)
        self.assertEqual(a.V1oc, 5 / s, "V1oc incorrect.")
        self.assertEqual(a.V2oc, 5 / s, "V2oc incorrect.")
        self.assertEqual(a.I1sc, 0.5 / s, "I1sc incorrect.")
        self.assertEqual(a.I2sc, -0.5 / s, "I2sc incorrect.")


    def test_LSection(self):
        """Mcircuits: check LSection

        """
        a = LSection(R(10), R(30))

        self.assertEqual(a.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(a.Igain12, -1, "Igain12 incorrect.")
        self.assertEqual(a.Z1oc, 40, "Z1oc incorrect.")
        self.assertEqual(a.Z2oc, 30, "Z2oc incorrect.")
        self.assertEqual(a.Z1sc, 10, "Z1sc incorrect.")
        self.assertEqual(a.Z2sc, 7.5, "Z2sc incorrect.")
        self.assertEqual(a.Vgain(1, 2), 0.75, "Vgain 1->2 incorrect.")
        self.assertEqual(a.Vgain(2, 1), 1, "Vgain 2->1 incorrect.")
        self.assertEqual(a.Igain(1, 2), -1, "Igain 1->2 incorrect.")
        self.assertEqual(a.Igain(2, 1), -0.75, "Igain 2->1 incorrect.")
        self.assertEqual(a.Vresponse(V(1), 1, 2), 0.75 / s, "Vresponse 1->2 incorrect.")
        self.assertEqual(a.Vresponse(V(1), 2, 1), 1 / s, "Vresponse 2->1 incorrect.")
        self.assertEqual(a.Iresponse(I(1), 1, 2), -1 / s, "Iresponse 1->2 incorrect.")
        self.assertEqual(a.Iresponse(I(1), 2, 1), -0.75 / s, "Iresponse 2->1 incorrect.")

        b = a.load(R(30))
        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z, 25, "R loaded Z incorrect.")
        self.assertEqual(b.V, 0, "R loaded V incorrect.")

        c = a.load(R(30) + V(5))
        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(c.Z, 25, "T loaded Z incorrect.")
        self.assertEqual(c.V, 2.5 / s, "T loaded V incorrect.")


        d = LSection(R(10), R(30) + V(6))
        self.assertEqual(d.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(d.Igain12, -1, "Igain12 incorrect.")
        self.assertEqual(d.Z1oc, 40, "Z1oc incorrect.")
        self.assertEqual(d.Z2oc, 30, "Z2oc incorrect.")
        self.assertEqual(d.V1oc, 6 / s, "V1oc incorrect.")
        self.assertEqual(d.V2oc, 6 / s, "V2oc incorrect.")
        self.assertEqual(d.I1sc, 0 / s, "I2sc incorrect.")
        self.assertEqual(d.I2sc, -0.2 / s, "I2sc incorrect.")
        

        e = d.chain(Series(R(10)))
        self.assertEqual(e.Z1oc, 40, "Z1oc incorrect.")
        self.assertEqual(e.Z2oc, 40, "Z2oc incorrect.")
        self.assertEqual(e.Vgain12, 0.75, "Vgain12 incorrect.")
        self.assertEqual(e.Igain12, -0.75, "Igain12 incorrect.")


    def test_Shunt_parallel(self):
        """Mcircuits: check Shunts in parallel

        """

        a = Shunt(R(10) + V(5))
        b = a.chain(Shunt(R(30)))
        self.assertEqual(b.Z1oc, 7.5, "Z1oc incorrect.")
        self.assertEqual(b.Z2oc, 7.5, "Z2oc incorrect.")
        self.assertEqual(b.V1oc, 3.75 / s, "V1oc incorrect.")
        self.assertEqual(b.V2oc, 3.75 / s, "V2oc incorrect.")


    def test_Series_series(self):
        """Mcircuits: check Series in series

        """

        a = Series(R(10) + V(5))
        b = a.chain(Series(R(30)))
        self.assertEqual(b.Y1sc, 0.025, "Y1sc incorrect.")
        self.assertEqual(b.Y2sc, 0.025, "Y2sc incorrect.")
        self.assertEqual(b.I1sc, -0.125 / s, "I1sc incorrect.")
        self.assertEqual(b.I2sc, 0.125 / s, "I2sc incorrect.")

        c = a.chain(a)



    def test_load(self):
        """Mcircuits: check load

        """

        a = Shunt(R(10) + V(5))
        b = a.load(R(30))

        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z, 7.5, "Shunt loaded R incorrect Z.")
        

    def test_open_circuit(self):
        """Mcircuits: check open_circuit

        """

        a = Shunt(R(10) + V(5))
        b = a.open_circuit()

        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z, 10, "incorrect Z.")
        self.assertEqual(b.V, 5 / s, "incorrect V.")


    def test_short_circuit(self):
        """Mcircuits: check short_circuit

        """

        a = Series(R(10) + V(5))
        b = a.short_circuit()

        self.assertEqual(type(b), Norton, "type incorrect.")
        self.assertEqual(b.Z, 10, "incorrect Z.")
        self.assertEqual(b.V, 5 / s, "incorrect V.")


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
        
