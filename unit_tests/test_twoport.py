from lcapy import *
import unittest
import sympy as sym
from msignal.mrf import MRF

s = sym.var('s')

class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy

    """

    def test_LSection(self):
        """Lcapy: check LSection

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
        self.assertEqual(a.Vresponse(Vdc(1), 1, 2), 0.75 / s, "Vresponse 1->2 incorrect.")
        self.assertEqual(a.Vresponse(Vdc(1), 2, 1), 1 / s, "Vresponse 2->1 incorrect.")
        self.assertEqual(a.Iresponse(Idc(1), 1, 2), -1 / s, "Iresponse 1->2 incorrect.")
        self.assertEqual(a.Iresponse(Idc(1), 2, 1), -0.75 / s, "Iresponse 2->1 incorrect.")

        b = a.load(R(30))
        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(b.Z, 25, "R loaded Z incorrect.")
        self.assertEqual(b.V, 0, "R loaded V incorrect.")

        c = a.load(R(30) + Vdc(5))
        self.assertEqual(type(b), Thevenin, "type incorrect.")
        self.assertEqual(c.Z, 25, "T loaded Z incorrect.")
        self.assertEqual(c.V, 2.5 / s, "T loaded V incorrect.")


        d = LSection(R(10), R(30) + Vdc(6))
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
        """Lcapy: check Shunts in parallel

        """

        a = Shunt(R(10) + Vdc(5))
        b = a.chain(Shunt(R(30)))
        self.assertEqual(b.Z1oc, 7.5, "Z1oc incorrect.")
        self.assertEqual(b.Z2oc, 7.5, "Z2oc incorrect.")
        self.assertEqual(b.V1oc, 3.75 / s, "V1oc incorrect.")
        self.assertEqual(b.V2oc, 3.75 / s, "V2oc incorrect.")


    def test_Series_series(self):
        """Lcapy: check Series in series

        """

        a = Series(R(10) + Vdc(5))
        b = a.chain(Series(R(30)))
        self.assertEqual(b.Y1sc, 0.025, "Y1sc incorrect.")
        self.assertEqual(b.Y2sc, 0.025, "Y2sc incorrect.")
        self.assertEqual(b.I1sc, -0.125 / s, "I1sc incorrect.")
        self.assertEqual(b.I2sc, 0.125 / s, "I2sc incorrect.")

        c = a.chain(a)


    def test_LSection_models(self):
        """Lcapy: check LSection models

        """
        a = LSection(R(10) + Vdc(5), R(20))

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
        
