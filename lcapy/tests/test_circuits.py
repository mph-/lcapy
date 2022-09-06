from lcapy import *
from lcapy.texpr import TimeDomainVoltage, TimeDomainCurrent
from lcapy.sexpr import (LaplaceDomainVoltage, LaplaceDomainTransferFunction,
                         LaplaceDomainImpedance, LaplaceDomainAdmittance)
import unittest
import sympy as sym
import os


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy circuits

    """

    def assertEqual2(self, ans1, ans2, comment):

        ans1 = ans1.canonical()
        ans2 = ans2.canonical()

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            ans1.pprint()
            ans2.pprint()
            raise AssertionError(e)

    def test_RC1(self):
        """Check RC network

        """
        a = Circuit()
        a.add('R1 1 2')
        a.add('C1 2 0')

        self.assertEqual2(a.impedance(1, 2), R(
            'R1').impedance, "Z incorrect for R1.")
        self.assertEqual2(a.impedance(2, 0), C(
            'C1').impedance, "Z incorrect for C1.")
        self.assertEqual2(
            a.impedance(1, 0), (R('R1') + C('C1')).impedance, "Z incorrect for R1 + C1.")

        self.assertEqual2(a.admittance(1, 2), R('R1').Y, "Y incorrect for R1.")
        self.assertEqual2(a.admittance(2, 0), C('C1').Y, "Y incorrect for C1.")
        self.assertEqual2(
            a.admittance(1, 0), (R('R1') + C('C1')).Y, "Y incorrect for R1 + C1.")
        self.assertEqual2(a.Isc(1, 0), I(0).Isc, "Isc incorrect")

        self.assertEqual(a.R1.Z, impedance('R1'), "Z incorrect")
        self.assertEqual(a.R1.R, impedance('R1'), "R incorrect")
        self.assertEqual(a.R1.X, 0, "X incorrect")
        self.assertEqual(a.C1.Y, admittance('s * C1'), "Y incorrect")
        self.assertEqual(a.C1.G, 0, "G incorrect")
        # Use ac to force jomega form.
        self.assertEqual(a.ac().C1.B, admittance('-omega * C1'), "B incorrect")

    def test_VRC1(self):
        """Check VRC circuit

        """
        a = Circuit()
        a.add('V1 1 0 {V1 / s}')
        a.add('R1 1 2')
        a.add('C1 2 0')

        # Note, V1 acts as a short-circuit for the impedance/admittance
        self.assertEqual2(
            a.impedance(1, 2), (R('R1') | C('C1')).Z, "Z incorrect across R1")
        self.assertEqual2(
            a.impedance(2, 0), (R('R1') | C('C1')).Z, "Z incorrect across C1")
        self.assertEqual2(a.impedance(1, 0), R(0).Z, "Z incorrect across V1")

        self.assertEqual2(
            a.admittance(1, 2), (R('R1') | C('C1')).Y, "Y incorrect across R1")
        self.assertEqual2(
            a.admittance(2, 0), (R('R1') | C('C1')).Y, "Y incorrect across C1")
        # This has a non-invertible A matrix.
        # self.assertEqual2(a.admittance(1, 0), R(0).Y, "Y incorrect across V1")

        self.assertEqual2(a.Voc(1, 0).s, V('V1 / s').Voc.s,
                          "Voc incorrect across V1")
        self.assertEqual(a.is_IVP, False, "Initial value problem incorrect")
        self.assertEqual(a.is_dc, False, "DC incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")

    def test_VRL1(self):
        """Check VRL circuit

        """
        a = Circuit()
        a.add('V1 1 0 {V1 / s}')
        a.add('R1 1 2')
        a.add('L1 2 0 L1 0')

        # This currently fails due to two symbols of the same name
        # having different assumptions.

        # Note, V1 acts as a short-circuit for the impedance/admittance
        self.assertEqual2(
            a.thevenin(1, 2).Voc, a.Voc(1, 2), "incorrect thevenin voltage")
        self.assertEqual2(
            a.thevenin(1, 2).Z, a.impedance(1, 2), "incorrect thevenin impedance")
        self.assertEqual2(
            a.norton(1, 2).Isc, a.Isc(1, 2), "incorrect norton current")
        self.assertEqual2(
            a.norton(1, 2).Y, a.admittance(1, 2), "incorrect norton admittance")
        self.assertEqual2(
            a.impedance(1, 2), (R('R1') | L('L1')).Z, "Z incorrect across R1")
        self.assertEqual2(
            a.impedance(2, 0), (R('R1') | L('L1')).Z, "Z incorrect across L1")
        self.assertEqual2(a.impedance(1, 0), R(0).Z, "Z incorrect across V1")

        self.assertEqual2(
            a.admittance(1, 2), (R('R1') | L('L1')).Y, "Y incorrect across R1")
        self.assertEqual2(
            a.admittance(2, 0), (R('R1') | L('L1')).Y, "Y incorrect across L1")
        # This has a non-invertible A matrix.
        # self.assertEqual2(a.admittance(1, 0), R(0).Y, "Y incorrect across V1")

        self.assertEqual2(a.Voc(1, 0), V('V1' / s).Voc,
                          "Voc incorrect across V1")
        self.assertEqual(a.is_IVP, True, "Initial value problem incorrect")
        self.assertEqual(a.is_dc, False, "DC incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")

    def test_IR1(self):
        """Check IR circuit

        """
        a = Circuit()
        a.add('I1 1 0 2')
        a.add('R1 1 0 1')

        self.assertEqual2(a.R1.V, V(2).Voc, "Incorrect voltage")

        self.assertEqual2(a[1].V, V(2).Voc, "Incorrect node voltage")

    def test_VCVS1(self):
        """Check VCVS

        """
        a = Circuit()
        a.add('V1 1 0 2')
        a.add('R1 1 0 1')
        a.add('E1 2 0 1 0 3')
        a.add('R2 2 0 1')

        self.assertEqual2(a.R2.V, V(6).Voc, "Incorrect voltage")

    def test_VCCS1(self):
        """Check VCCS

        """
        a = Circuit()
        a.add('V1 1 0 2')
        a.add('R1 1 0 1')
        a.add('G1 2 0 1 0 3')
        a.add('R2 2 0 1')

        self.assertEqual2(a.R2.V, V(6).Voc, "Incorrect voltage")

    def test_CCCS1(self):
        """Check CCCS

        """
        a = Circuit()
        a.add('V1 1 0 10')
        a.add('R1 1 2 2')
        a.add('V2 2 0 0')
        a.add('F1 3 0 V2 2')
        a.add('R2 3 0 1')

        self.assertEqual2(a.R2.V, V(10).Voc, "Incorrect voltage")

    def test_CCVS1(self):
        """Check CCVS

        """
        a = Circuit()
        a.add('V1 1 0 10')
        a.add('R1 1 2 2')
        a.add('V2 2 0 0')
        a.add('H1 3 0 V2 2')
        a.add('R2 3 0 1')

        self.assertEqual2(a.R2.V, V(10).Voc, "Incorrect voltage")

    def test_CCVS2(self):
        """Check CCVS

        """
        a = Circuit()
        a.add('H1 1 0 H1 -4')
        t = a.thevenin(1, 0)
        self.assertEqual(t.V, 0, "Incorrect Thevenin voltage")
        self.assertEqual(t.Z, impedance(-4), "Incorrect Thevenin impedance")
        self.assertEqual(a.H1.Voc, 0, "Incorrect cpt voltage")
        self.assertEqual(a.H1.Z, 0, "Incorrect cpt impedance")

    def test_V1(self):
        """Test V1"""

        a = Circuit()
        a.add('V1 1 0 10')

        self.assertEqual(a.V1.V.dc, voltage(10), "Incorrect voltage")

    def test_VRL1_dc(self):
        """Check VRL circuit at dc

        """

        a = Circuit()
        a.add('V1 1 0 dc')
        a.add('R1 1 2')
        a.add('L1 2 0')
        self.assertEqual(a.is_IVP, False, "Initial value problem incorrect")
        self.assertEqual(a.is_dc, True, "DC incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")

    def test_VRL1_dc2(self):
        """Check VRL circuit at dc but with initial conditions

        """

        # TODO: This currently fails due to two symbols of the same name
        # having different assumptions.

        a = Circuit()
        a.add('V1 1 0 {V1 + 1}')
        a.add('R1 1 2')
        a.add('L1 2 0 L1 {(V1 + 1) / R1}')
        # This tests if symbols are converted to the defined ones.
        self.assertEqual2(a.L1.v, V(0).Voc.s.inverse_laplace(),
                          "Incorrect time domain voltage")
        v = LaplaceDomainVoltage('(V1+1)/s', dc=False).inverse_laplace()
        self.assertEqual2(a.R1.v, v,
                          "Incorrect time domain voltage")
        self.assertEqual(a.is_IVP, True, "Initial value problem incorrect")
        self.assertEqual(a.is_dc, False, "DC incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")

    def test_VRL1_dc3(self):
        """Check VRL circuit at dc but with initial conditions

        """

        # TODO: This currently fails due to two symbols of the same name
        # having different assumptions.

        a = Circuit()
        a.add('V1 1 0 V1')
        a.add('R1 1 2')
        a.add('L1 2 0 L1 0')
        self.assertEqual(a.is_dc, False, "DC incorrect")
        self.assertEqual(a.is_ac, False, "AC incorrect")

    def test_VRL1_ac(self):
        """Check VRL circuit at ac

        """

        a = Circuit()
        a.add('V1 1 0 ac')
        a.add('R1 1 2')
        a.add('L1 2 0')
        self.assertEqual(a.is_IVP, False, "Initial value problem incorrect")
        self.assertEqual(a.is_dc, False, "DC incorrect")
        self.assertEqual(a.is_ac, True, "AC incorrect")
        self.assertEqual(a.R1.I, a.L1.I, "currents different")
        self.assertEqual(a.V1.I, a.L1.I, "currents different")

    def test_VRC2(self):
        """Check VRC circuit with arbitrary s-domain source

        """

        a = Circuit()
        a.add('V1 1 0 {V(s)}')
        a.add('R1 1 2')
        a.add('C1 2 0 C1 0')
        H = a[2].V.s / a[1].V.s

        self.assertEqual2(H, 1 / (s * 'R1' * 'C1' + 1),  "Incorrect ratio")

    def test_causal1(self):
        """Test VRL circuit causality"""

        a = Circuit()
        a.add('V1 1 0 {4 + 2 * u(t)}; down')
        a.add('R1 1 2 2; right=2')
        a.add('L1 2 3 2; down')
        a.add('W 0 3; right')

        self.assertEqual(a.sub['s'].is_causal, True, "Causal incorrect")
        self.assertEqual2(a.L1.v, voltage(
            2 * exp(-t) * u(t)), "L current incorrect")

    def test_VR1_ac2(self):
        """Check VR circuit at ac for angular frequency 1

        """

        # This can be solved in the time-domain
        a = Circuit()
        a.add('V1 1 0 ac 5 0 1')
        a.add('R1 1 0 1')
        self.assertEqual(a.is_IVP, False, "Initial value problem incorrect")
        self.assertEqual(a.is_dc, False, "DC incorrect")
        self.assertEqual(a.is_ac, True, "AC incorrect")
        self.assertEqual(a.is_causal, False, "Causal incorrect")
        self.assertEqual(a.is_time_domain, True, "Time domain incorrect")
        self.assertEqual(a.V1.v, TimeDomainVoltage(
            '5*cos(t)'), "V1 voltage incorrect")
        self.assertEqual(a.R1.v, TimeDomainVoltage(
            '5*cos(t)'), "R1 voltage incorrect")
        self.assertEqual(a.V1.i, TimeDomainCurrent(
            '5*cos(t)'), "V1 current incorrect")
        self.assertEqual(a.R1.i, TimeDomainCurrent(
            '5*cos(t)'), "R1 current incorrect")

    def test_VRL1_ac2(self):
        """Check VRL circuit at ac for angular frequency 1

        """

        a = Circuit()
        a.add('V1 1 0 ac 5 0 1')
        a.add('R1 1 2 3')
        a.add('L1 2 0 4')
        self.assertEqual(a.is_IVP, False, "Initial value problem incorrect")
        self.assertEqual(a.is_dc, False, "DC incorrect")
        self.assertEqual(a.is_ac, True, "AC incorrect")
        self.assertEqual(a.is_causal, False, "Causal incorrect")
        self.assertEqual(a.is_time_domain, False, "Time domain incorrect")
        self.assertEqual(a.V1.v, TimeDomainVoltage(
            '5*cos(t)'), "V1 voltage incorrect")
        self.assertEqual(a.R1.i, TimeDomainCurrent(
            '(4*sin(t)+3*cos(t))/5'), "R1 current incorrect")

    def test_VRC_ivp(self):
        """Check VRC IVP"""

        a = Circuit("""
        V 1 0 dc; down
        R 1 2; right
        C 2 3 C V0; down
        W 0 3; right""")

        self.assertEqual(a.is_IVP, True, "Initial value problem incorrect")
        self.assertEqual(a.R.I, a.C.I, "R + C current different")
        self.assertEqual(a.V.I, a.C.I, "V + C current different")
        self.assertEqual(a.V.V,  a.R.V + a.C.V, "KVL fail")

        a = Circuit("""
        V 1 0 dc; down
        R 2 3; right
        C 1 2 C V0; down
        W 0 3; right""")

        self.assertEqual(a.is_IVP, True, "Initial value problem incorrect")
        self.assertEqual(a.R.I, a.C.I, "R + C current different")
        self.assertEqual(a.V.I, a.C.I, "V + C current different")
        self.assertEqual(a.V.V,  a.R.V + a.C.V, "KVL fail")

    def test_VRL_ivp(self):
        """Check VRL IVP"""

        a = Circuit("""
        V 1 0 dc; down
        R 1 2; right
        L 2 3 L I0; down
        W 0 3; right""")

        self.assertEqual(a.is_IVP, True, "Initial value problem incorrect")
        self.assertEqual(a.R.I, a.L.I, "R + L current different")
        self.assertEqual(a.V.I, a.L.I, "V + L current different")
        self.assertEqual(a.V.V,  a.R.V + a.L.V, "KVL fail")

        a = Circuit("""
        V 1 0 dc; down
        R 2 3; right
        L 1 2 L I0; down
        W 0 3; right""")

        self.assertEqual(a.is_IVP, True, "Initial value problem incorrect")
        self.assertEqual(a.R.I, a.L.I, "R + L current different")
        self.assertEqual(a.V.I, a.L.I, "V + L current different")
        self.assertEqual(a.V.V,  a.R.V + a.L.V, "KVL fail")

    def test_RL_ivp(self):
        """Check RL IVP"""

        a = Circuit("""
        L 1 0 L i0
        R 1 0""")

        self.assertEqual(a.is_IVP, True, "Initial value problem incorrect")
        # Note, a positive current through an inductor is from the
        # positive node to the negative node.
        self.assertEqual(a.R.I, -a.L.I, "R + L current different")
        self.assertEqual(a.R.V, a.L.V, "R + L voltage different")
        self.assertEqual(a.L.I(s), current(
            'i0 / (s + R / L)'), "L current wrong")

    def test_RC_ivp(self):
        """Check RC IVP"""

        a = Circuit("""
        C 1 0 C v0
        R 1 0""")

        self.assertEqual(a.is_IVP, True, "Initial value problem incorrect")
        self.assertEqual(a.R.I, -a.C.I, "R + C current different")
        self.assertEqual(a.R.V, a.C.V, "R + C voltage different")
        self.assertEqual(a.C.I(s), current(
            '-v0 / (s * R + 1 / C)'), "C current wrong")

    def test_sub(self):
        """Test sub-circuits"""

        a = Circuit("""
        V1 1 0 {u(t)}
        R1 1 2
        L1 2 0""")

        self.assertEqual(a.ac().R1.V, 0, "AC model incorrect")
        self.assertEqual(a.dc().R1.V, 0, "DC model incorrect")

    def test_IV_series(self):
        """Test IV series"""

        a = Circuit("""
        V 2 1 dc
        I 1 0 dc
        R 2 0""")

        self.assertEqual(a.R.v, voltage('I * R'), "R voltage incorrect")
        self.assertEqual(a.R.i, current('I'), "R current incorrect")

    def test_IV_parallel(self):
        """Test IV parallel"""

        a = Circuit("""
        V 1 0 dc
        I 1 0 dc
        R 1 0""")

        self.assertEqual(a.R.v, voltage('V'), "R voltage incorrect")
        self.assertEqual(a.R.i, current('V / R'), "R current incorrect")

    def test_oneport(self):
        """Test oneport"""

        a = Circuit("""
        V 2 0 s
        L 2 1
        R 1 0""")

        th = a[1].thevenin()
        no = a[1].norton()

        self.assertEqual(th.Voc, no.Voc, "Voc incorrect")
        self.assertEqual(th.Isc, no.Isc, "Isc incorrect")
        self.assertEqual(th.Z, no.Z, "Z incorrect")
        self.assertEqual(th.Y, no.Y, "Y incorrect")

        th = a[2].thevenin(1)
        no = a[2].norton(1)

        self.assertEqual(th.Voc, no.Voc, "Voc incorrect")
        self.assertEqual(th.Isc, no.Isc, "Isc incorrect")
        self.assertEqual(th.Z, no.Z, "Z incorrect")
        self.assertEqual(th.Y, no.Y, "Y incorrect")

        th = a.L.thevenin()
        no = a.L.norton()

        self.assertEqual(th.Voc, no.Voc, "Voc incorrect")
        self.assertEqual(th.Isc, no.Isc, "Isc incorrect")
        self.assertEqual(th.Z, no.Z, "Z incorrect")
        self.assertEqual(th.Y, no.Y, "Y incorrect")

        th = a.thevenin(2, 1)
        no = a.norton(2, 1)

        self.assertEqual(th.Voc, no.Voc, "Voc incorrect")
        self.assertEqual(th.Isc, no.Isc, "Isc incorrect")
        self.assertEqual(th.Z, no.Z, "Z incorrect")
        self.assertEqual(th.Y, no.Y, "Y incorrect")

    def test_directive(self):
        """Test directive"""

        # This has a deliberate empty line.
        a = Circuit("""
        V 2 0 10

        R 2 0 1""")

        self.assertEqual(a.R.i, current(10), "i incorrect")

    def test_K(self):
        """Test mutual inductances"""

        a = Circuit("""
        V1 1 0 ac; down
        W 1 1_1; right
        L1 1_1 0_1 L; down
        W 0 0_1; right
        W 0_1 0_2; right
        W 0_1 0_2; right
        L2 2 0_2 L; down
        W 2 2_1; right
        R1 2_1 0_3; down
        W 0_2 0_3; right
        K L1 L2 1; invisible""")

        self.assertEqual(a.R1.v, voltage('V1 * cos(omega_0 * t) * abs(omega_0) / omega_0'),
                         "coupling incorrect")

    def test_kill(self):
        """Test kill"""

        a = Circuit("""
        V1 1 2 6
        V2 2 0 4
        R 1 0 2""")

        self.assertEqual(a.R.V, voltage(10), "incorrect series voltage")

        b = a.kill('V1')

        self.assertEqual(b.R.V, voltage(4), "incorrect voltage with V1 killed")

        c = a.kill('V2')

        self.assertEqual(c.R.V, voltage(6), "incorrect voltage with V2 killed")

        d = a.kill_except('V1')

        self.assertEqual(d.R.V, voltage(
            6), "incorrect voltage with all killed except V1")

    def test_TF(self):
        """Test TF"""

        a = Circuit("""
        TF 2 0 1 0 k
        R 2 0""")

        self.assertEqual(a.impedance(1, 0), impedance(
            'R / k**2'), "incorrect impedance")

    def test_TVtriode(self):
        """ Test thermionic valve triode"""

        if False:
            a = Circuit("""
            TVtriode 0 1 2 100 1.6e-3 3e-12 3e-12; up
            """).draw("triode.png")
            self.assertTrue(os.path.exists('triode.png'))
            os.remove('triode.png')

    def test_params(self):
        """Test params"""

        a = Circuit("""
        R1 1 2;
        R2 2 0;
        R3 2 3""")

        A = a.Aparams(3, 0, 1, 0)

        # TODO, FIXME without passing matrix elements through expr
        self.assertEqual(expr(A[0, 0]), expr('(R2 + R3) / R2'), "A11")
        self.assertEqual(expr(A[0, 1]), expr('R1 + R1 * R3 / R2 + R3'), "A12")
        self.assertEqual(expr(A[1, 0]), expr('1 / R2'), "A21")
        self.assertEqual(expr(A[1, 1]), expr('(R1 + R2) / R2'), "A22")

        Z = a.Zparams(3, 0, 1, 0)

        self.assertEqual(expr(Z[0, 0]), expr('R2 + R3'), "Z11")
        self.assertEqual(expr(Z[0, 1]), expr('R2'), "Z12")
        self.assertEqual(expr(Z[1, 0]), expr('R2'), "Z21")
        self.assertEqual(expr(Z[1, 1]), expr('R1 + R2'), "Z22")

    def test_in_series(self):
        """Test in_series"""

        a = Circuit("""
        V 1 0 dc
        R1 1 2
        R2 2 3
        R3 3 0
        R4 3 4
        R5 4 0""")

        self.assertEqual(set(a.in_series('R1')),
                         set(('V', 'R2')), "in_series(R1)")
        self.assertEqual(set(a.in_series('R3')), set(), "in_series(R3)")
        self.assertEqual(set(a.in_series('R4')), set(('R5',)), "in_series(R4)")

    def test_in_parallel(self):
        """Test in_parallel"""

        a = Circuit("""
        V 1 0 dc
        R1 1 0
        R2 1 0
        R3 1 2
        R4 2 0
        R5 2 0""")

        self.assertEqual(set(a.in_parallel('R1')),
                         set(('V', 'R2')), "in_parallel(R1)")
        self.assertEqual(set(a.in_parallel('R3')), set(), "in_parallel(R3)")
        self.assertEqual(set(a.in_parallel('R4')),
                         set(('R5', )), "in_parallel(R4)")

    def test_simplify(self):
        """Test simplify"""

        a = Circuit("""
        R1 1 2
        R2 2 3""")

        b = a.simplify()

        self.assertEqual(b.impedance(
            1, 3), a.impedance(1, 3), "simplify series")

        a = Circuit("""
        R1 1 2
        R2 1 2""")

        b = a.simplify()

        self.assertEqual(b.impedance(1, 2), a.impedance(
            1, 2), "simplify parallel")

    def test_VRL1_super(self):
        """Check VRL circuit

        """
        a = Circuit()
        a.add('V1 1 0 {2 * u(t)}; down')
        a.add('V2 1 4 4; right')
        a.add('R1 4 2 2; right=2')
        a.add('L1 2 3 2; down')
        a.add('W 0 3; right')

    def test_super_nodes(self):
        """Test super_nodes"""

        a = Circuit("""
        V1 1 2;
        R1 2 3:
        R2 1 3""")

        self.assertEqual(a.super_nodes, [['1', '2']], "super_nodes")

    def test_transfer(self):
        """Test transfer"""

        a = Circuit("""
        V 2 0
        L 2 1
        R 1 0""")

        H1 = a.L.transfer('R')
        H2 = a.transfer(2, 1, 1, 0)

        self.assertEqual(H1, H2, "transfer cpt/nodes")
        self.assertEqual(H1, -1, "transfer")

    def test_replace(self):
        """Test replace"""

        a = Circuit("""
        V 2 0
        L 2 1
        R 1 0""")

        b = a.replace('V', 'W')
        c = a.replace('V', 'W 2 0')

        self.assertEqual(a.has('V'), True, "has V")
        self.assertEqual(b.has('V'), False, "has V")
        self.assertEqual(c.has('V'), False, "has V")

    def test_connected(self):
        """Test connected"""

        a = Circuit("""
        V 2 0
        L 2 1
        R 1 0""")

        self.assertEqual(a.is_connected, True, "is_connected")

        a.add("R2 3 4")

        self.assertEqual(a.is_connected, False, "not is_connected")

    def test_netlist_subs(self):
        """Test netlist subs"""

        a = Circuit()
        a.add('V1 1 0; down')

        b = a.subs({'V1': f})

        self.assertEqual(b.V1.Voc, expr(f), "netlist subs")

    def test_cpt(self):
        """Test cpt"""

        a = Circuit()
        a.add('V1 1 0')
        a.add('E1 2 0 opamp 1 0')
        a.add('R1 2 3')

        self.assertEqual(a.V1.is_source, True, "V is_source")
        self.assertEqual(a.V1.is_dependent_source,
                         False, "V is_dependent_source")
        self.assertEqual(a.V1.is_independent_source,
                         True, "V is_independent_source")
        self.assertEqual(a.V1.is_resistor, False, "V is_resistor")

        self.assertEqual(a.E1.is_source, True, "E is_source")
        self.assertEqual(a.E1.is_dependent_source,
                         True, "E is_dependent_source")
        self.assertEqual(a.E1.is_independent_source,
                         False, "E is_independent_source")
        self.assertEqual(a.E1.is_resistor, False, "E is_resistor")

        self.assertEqual(a.R1.is_source, False, "E is_source")
        self.assertEqual(a.R1.is_dependent_source,
                         False, "R is_dependent_source")
        self.assertEqual(a.R1.is_independent_source,
                         False, "R is_independent_source")
        self.assertEqual(a.R1.is_resistor, True, "E is_resistor")

    def test_TP(self):
        """Test TP"""

        n = Series(R(5)).chain(Shunt(R(10)))
        b = n.Bmodel
        c = b.circuit()
        c.add('V 1 0 15; down')
        self.assertEqual(c[3].V, 10, "Bmodel")

        d = n.Ymodel
        e = d.circuit()
        e.add('V 1 0 15; down')
        self.assertEqual(e[3].V, 10, "Ymodel")

    def test_transfer(self):
        """Check transfer function

        """

        a = Circuit()
        a.add('R1 1 0 1')
        a.add('R2 1 2 2')
        a.add('C1 2 0 1')

        H = a.transfer(1, 0, 2, 0)
        self.assertEqual2(H, 1 / (2 * s + 1), "Incorrect transfer function")
        h = H.inverse_laplace()
        self.assertEqual2(h, exp(-t / 2) * Heaviside(t) / 2,
                          "Incorrect impulse response")
        H1 = a.transfer((1, 0), (2, 0))
        self.assertEqual2(H1, H, "Incorrect transfer function")
        H2 = a.transfer('R1', 'C1')
        self.assertEqual2(H2, H, "Incorrect transfer function")

    def test_transfer(self):
        """Check transfer functions

        """

        a = Circuit("""
        P1 1 0
        R1 1 2
        R2 2 0
        R3 2 3
        P2 3 0""")

        Av = a.voltage_gain('P1', 'P2')
        self.assertEqual2(Av, LaplaceDomainTransferFunction('R2 / (R1 + R2)'),
                          'voltage_gain')

        Ai = a.current_gain('P1', 'P2')
        self.assertEqual2(Ai, LaplaceDomainTransferFunction('-R2 / (R2 + R3)'),
                          'current_gain')

        Z = a.transimpedance('P1', 'P2')
        self.assertEqual2(Z, LaplaceDomainImpedance('R2'),
                          'transimpedance')

        Y = a.transadmittance('P1', 'P2')
        self.assertEqual2(Y, LaplaceDomainAdmittance('-R2 / (R1 * R2 + R1 * R3 + R2 * R3)'),
                          'transadmittance')

    def test_expand(self):
        """Test expand"""

        a = Circuit('E1 3 0 opamp 1 2\n')
        b = a.expand()

        self.assertEqual(str(b), 'E__E1 3 0 1 2 E1 0', 'opamp expand')

    def test_R_simplify(self):

        a = Circuit("""
        R1 1 2 10
        R2 2 3 5
        R3 3 4 15""")

        self.assertEqual(a.simplify().Rt1.R, 30, 'series R')
        self.assertEqual(a.simplify_series().Rt1.R, 30, 'series R')
        self.assertEqual(a.simplify_series(
            ignore=['R3']).Rt1.R, 15, 'series R')

        self.assertEqual(len(a.simplify().remove_dangling_wires().elements),
                         1, 'dangling wires')

        b = Circuit("""
        R1 1 2 20
        R2 1 2 60
        R3 1 2 15""")

        self.assertEqual(b.simplify().Rt1.R, 7.5, 'parallel R')
        self.assertEqual(b.simplify_parallel().Rt1.R, 7.5, 'parallel R')
        self.assertEqual(b.simplify_parallel(
            ignore=['R1']).Rt1.R, 12, 'parallel R')
        self.assertEqual(b.simplify_parallel(
            select=['R2', 'R3']).Rt1.R, 12, 'parallel R')

        c = Circuit("""
        R1 1 2 20
        R2 1 2 60
        R3 3 4 15""")

        self.assertEqual(c.simplify().Rt1.R, 15,
                         'parallel R with disconnected')

        self.assertEqual(len(c.remove_disconnected().elements),
                         2, 'disconnected')
