from lcapy import *
from lcapy.cexpr import ConstantDomainExpression, ConstantFrequencyResponseDomainExpression, ConstantTimeDomainExpression
from lcapy.fexpr import FourierDomainExpression
from lcapy.omegaexpr import AngularFourierDomainExpression
from lcapy.texpr import TimeDomainExpression, TimeDomainVoltage
from lcapy.sexpr import LaplaceDomainExpression, LaplaceDomainVoltage
import unittest


class LcapyCoreTester(unittest.TestCase):
    """Unit tests for lcapy core"""

    def assertEqual2(self, ans1, ans2, comment):

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            pprint(ans1)
            pprint(ans2)
            raise AssertionError(e)

    def test_Expr1(self):
        """Check Expr

        """
        a = expr(3)
        self.assertEqual2(a.N, expr(3), "N incorrect.")
        self.assertEqual2(a.D, expr(1), "D incorrect.")
        self.assertEqual2(a.real, expr(3), "real incorrect.")
        self.assertEqual2(a.imag, expr(0), "imag incorrect.")
        self.assertEqual2(a.magnitude, expr(3), "magnitude incorrect.")
        self.assertEqual2(a.phase, rad(0), "phase incorrect.")
        self.assertEqual2(a.phase_degrees, deg(0), "phase incorrect.")
        self.assertEqual2(a.sign, expr(1), "sign incorrect.")
        self.assertEqual2(-a.sign, expr(-1), "sign incorrect.")

        a = expr(-3)
        self.assertEqual2(a.N, expr(-3), "N incorrect.")
        self.assertEqual2(a.D, expr(1), "D incorrect.")
        self.assertEqual2(a.real, expr(-3), "real incorrect.")
        self.assertEqual2(a.imag, expr(0), "imag incorrect.")
        self.assertEqual2(a.magnitude, expr(3), "magnitude incorrect.")
        self.assertEqual2(a.phase, rad(pi), "phase incorrect.")
        self.assertEqual2(a.phase_degrees, deg(180), "phase incorrect.")
        self.assertEqual2(a.sign, expr(-1), "sign incorrect.")
        self.assertEqual2(-a.sign, expr(1), "sign incorrect.")

        a = expr(3j)
        self.assertEqual2(a.N, expr(3j), "N incorrect.")
        self.assertEqual2(a.D, expr(1), "D incorrect.")
        self.assertEqual2(a.real, expr(0), "real incorrect.")
        self.assertEqual2(a.imag, expr(3), "imag incorrect.")
        self.assertEqual2(a.magnitude, expr(3), "magnitude incorrect.")
        self.assertEqual2(abs(a), expr(3), "abs incorrect.")
        self.assertEqual2(a.phase, rad(pi / 2), "phase incorrect.")
        self.assertEqual2(a.phase_degrees, deg(90), "phase incorrect.")
        self.assertEqual2(a.sign, expr(j), "sign incorrect.")
        self.assertEqual2(-a.sign, expr(-j), "sign incorrect.")

        n = s + 2
        d = (s + 3) * (s + 4)
        e = exp(-5 * s)
        a = n * e / d
        self.assertEqual2(a.N, n * e, "N incorrect.")
        self.assertEqual2(a.D, d, "D incorrect.")

    def test_cExpr1(self):
        """Check cExpr1

        """
        a = cexpr('1')
        self.assertEqual2(a.evaluate(), 1, "evaluate incorrect.")
        self.assertEqual2(type(a), ConstantDomainExpression, "type incorrect.")

        Z = impedance(a)
        NZ = Z.N
        self.assertEqual2(
            type(NZ), ConstantFrequencyResponseDomainExpression, "N type incorrect for Z.")

        V = voltage(a)
        NV = V.N
        self.assertEqual2(type(NV), ConstantTimeDomainExpression,
                          "N type incorrect for V.")

    def test_sExpr1(self):
        """Check sExpr1

        """
        a = LaplaceDomainExpression('(s+2)/(s-2)')
        self.assertEqual2(a.N, LaplaceDomainExpression('s+2'), "N incorrect.")
        self.assertEqual2(a.D, LaplaceDomainExpression('s-2'), "D incorrect.")
        self.assertEqual2(a.Ndegree, 1, "N degree incorrect.")
        self.assertEqual2(a.Ddegree, 1, "D degree incorrect.")
        self.assertEqual2(a.degree, 1, "Degree incorrect.")
        self.assertEqual2(a.is_rational_function, True,
                          "is_rational_function incorrect.")
        self.assertEqual2(a.is_strictly_proper, False,
                          "is_strictly_proper incorrect.")

        self.assertEqual2(sorted(a.poles()), [2], "poles incorrect.")
        self.assertEqual2(sorted(a.zeros()), [-2], "zeros incorrect.")

        self.assertEqual2(a.partfrac(), 1 + 4 / (s - 2), "partfrac incorrect.")

        self.assertEqual(a.evaluate(1), -3.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate(-2), 0.0, "scalar evaluate incorrect.")

        self.assertEqual(a.inverse_laplace(causal=True), 4 * exp(2 * t)
                         * H(t) + DiracDelta(t), "inverse Laplace incorrect.")

        aw = a.subs(j * omega, causal=True)
        self.assertEqual2(aw.real, (omega**2 - 4) /
                          (omega**2 + 4), "real part incorrect.")
        self.assertEqual2(aw.imag, -4 * omega /
                          (omega**2 + 4), "imag part incorrect.")

    def test_sExpr2(self):
        """Check sExpr2

        """
        a = (s + 2) * (s + 3) / (s - 2)
        self.assertEqual2(a.N, (s + 2) * (s + 3), "N incorrect.")
        self.assertEqual2(a.D, s - 2, "D incorrect.")

        self.assertEqual2(sorted(a.poles()), [2], "poles incorrect.")
        self.assertEqual2(sorted(a.zeros()), sorted(
            [-3.0, -2.0]), "zeros incorrect.")

        self.assertEqual2(
            a.partfrac(), s + 7 + 20 / (s - 2), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), s + 7 + 20 / (s - 2), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")
        self.assertEqual2(
            a.canonical(), (s**2 + 5 * s + 6) / (s - 2), "canonical incorrect.")

        self.assertEqual(a.inverse_laplace(causal=True), 20 * exp(2 * t) * H(t) +
                         7 * DiracDelta(t) + DiracDelta(t, 1), "inverse Laplace incorrect.")

    def test_sExpr3(self):
        """Check sExpr3

        """
        a = (s**2 + 5 * s + 6) / (s - 2)
        self.assertEqual2(a.N, s ** 2 + 5 * s + 6, "N incorrect.")
        self.assertEqual2(a.D, s - 2, "D incorrect.")
        self.assertEqual2(a.Ndegree, 2, "N degree incorrect.")
        self.assertEqual2(a.Ddegree, 1, "D degree incorrect.")
        self.assertEqual2(a.degree, 2, "Degree incorrect.")
        self.assertEqual2(a.is_strictly_proper, False,
                          "is_strictly_proper incorrect.")

        self.assertEqual2(sorted(a.poles()), [2], "poles incorrect.")
        self.assertEqual2(sorted(a.zeros()), sorted(
            [-3, -2]), "zeros incorrect.")

        self.assertEqual2(
            a.partfrac(), s + 7 + 20 / (s - 2), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), s + 7 + 20 / (s - 2), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")
        self.assertEqual2(
            a.canonical(), (s**2 + 5 * s + 6) / (s - 2), "canonical incorrect.")
        self.assertEqual2(
            a.ZPK(), (s + 2) * (s + 3) / (s - 2), "ZPK incorrect.")

        self.assertEqual(a.inverse_laplace(causal=True), 20 * exp(2 * t) * H(t) +
                         7 * DiracDelta(t) + DiracDelta(t, 1), "inverse Laplace incorrect.")

    def test_sExpr4(self):
        """Check sExpr4

        """
        a = 1 / ((s - j) * (s + j))
        self.assertEqual2(a.N, 1, "N incorrect.")
        self.assertEqual2(a.D, (s - j) * (s + j), "D incorrect.")

        self.assertEqual2(set(a.poles()), set(
            expr([-j, j])), "poles incorrect.")
        self.assertEqual2(a.zeros(), {}, "zeros incorrect.")

        # This depends on if 2 * (s + j) is expanded to 2 * s + 2 * j
        # self.assertEqual2(
        # a.partfrac(), j / (2 * (s + j)) - j / (2 * (s - j)), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), 1 / (s**2 + 1), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), 1 / (s**2 + 1), "general incorrect.")
        self.assertEqual2(
            a.canonical(), 1 / (s**2 + 1), "canonical incorrect.")
        self.assertEqual2(
            a.ZPK(), 1 / ((s - j) * (s + j)), "ZPK incorrect.")

        self.assertEqual(a.inverse_laplace(causal=True), sin(t)
                         * H(t), "inverse Laplace incorrect.")

    def test_sExpr5(self):
        """Check sExpr5

        """
        a = 1 / ((s - 1j) * (s + 1j))
        self.assertEqual2(a.N, 1, "N incorrect.")
        self.assertEqual2(a.D, (s - j) * (s + j), "D incorrect.")

        self.assertEqual2(set(a.poles()), set(
            expr([-j, j])), "poles incorrect.")
        self.assertEqual2(len(a.zeros()), 0, "zeros incorrect.")

        # This depends on if 2 * (s + j) is expanded to 2 * s + 2 * j
        # self.assertEqual2(
        #    a.partfrac(), j / (2 * (s + j)) - j / (2 * (s - j)), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), 1 / (s**2 + 1), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), 1 / (s**2 + 1), "general incorrect.")
        self.assertEqual2(
            a.canonical(), 1 / (s**2 + 1), "canonical incorrect.")
        self.assertEqual2(
            a.ZPK(), 1 / ((s - j) * (s + j)), "ZPK incorrect.")

        self.assertEqual(a.inverse_laplace(causal=True), sin(t)
                         * H(t), "inverse Laplace incorrect.")

    def test_sExpr6(self):
        """Check sExpr6 (repeated poles)

        """
        a = 1 / ((s + 4) * (s + 4))
        self.assertEqual2(a.N, 1, "N incorrect.")
        self.assertEqual2(a.D, (s + 4)**2, "D incorrect.")

        self.assertEqual2(a.poles(), {-4: 2}, "poles incorrect.")
        self.assertEqual2(a.zeros(), {}, "zeros incorrect.")

        self.assertEqual2(
            a.mixedfrac(), 1 / (s**2 + 8 * s + 16), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), 1 / (s**2 + 8 * s + 16), "general incorrect.")
        self.assertEqual2(
            a.canonical(), 1 / (s**2 + 8 * s + 16), "canonical incorrect.")
        self.assertEqual2(
            a.ZPK(), 1 / ((s + 4)**2), "ZPK incorrect.")
        self.assertEqual(a.inverse_laplace(causal=True), t *
                         exp(-4 * t) * Heaviside(t), "inverse Laplace incorrect.")

    def test_sExpr7(self):
        """Check sExpr7 (delay)

        """
        a = LaplaceDomainExpression('(s+1)*exp(-3*s)/((s+3)*(s+4))')
        self.assertEqual(a.inverse_laplace(causal=True), (3 * exp(-4 * t + 12) -
                         2 * exp(-3 * t + 9)) * H(t - 3), "inverse Laplace incorrect.")

    def test_sExpr8(self):
        """Check sExpr8 (jomega)

        """
        a = expr('a')
        H = a * s / (a * s + 1)
        Hw = H(j * omega, causal=True)

        self.assertEqual(Hw.phase, rad(atan2(1, a * omega)), 'phase')

    def test_sExpr9(self):
        """Check sExpr9 (non-monic factor)

        """
        H = 1 / ((s + 1) * (3 * s + 2))
        # Note, one is used to convert 2 to a SymPy Integer and thus
        # 2 / 3 is not converted into a floating point value
        # before becoming a rational number.
        H2 = 1 / (s + one * 2 / 3) - 1 / (s + 1)
        self.assertEqual2(H.partfrac(), H2, "partfrac incorrect.")

    def test_wExpr1(self):
        """Check wExpr1

        """

        A = (j * omega + 3) / (j * omega - 4)
        self.assertEqual2(A.N, j * omega + 3, "N incorrect.")
        self.assertEqual2(A.D, j * omega - 4, "D incorrect.")
        self.assertEqual2(A.real, (omega**2 - 12) /
                          (omega**2 + 16), "real incorrect.")
        self.assertEqual2(A.imag, -7 * omega /
                          (omega**2 + 16), "imag incorrect.")

    def test_tExpr1(self):
        """Check tExpr1

        """
        a = t**2

        self.assertEqual(a.evaluate(2), 4.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate((2, 3))[
                         1], 9.0, "vector evaluate incorrect.")

    def test_expr(self):
        """Check expr
        """

        e = expr('rect(t)')
        self.assertEqual(e, rect(t), "expr('rect(t)'")
        e = expr('tri(t)')
        self.assertEqual(e, tri(t), "expr('tri(t)'")

    def test_step(self):
        """Check step

        """
        a = u(t)

        self.assertEqual(a.evaluate(2), 1.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate((-2, 2))
                         [1], 1.0, "vector evaluate incorrect.")
        self.assertEqual(a.evaluate((-2, 2))
                         [0], 0.0, "vector evaluate incorrect.")
        self.assertEqual(a.laplace(), 1 / s, "Laplace transform incorrect.")

    def test_delta(self):
        """Check delta

        """
        a = delta(t)

        self.assertEqual(a.evaluate(2), 0, "scalar evaluate incorrect.")
        self.assertEqual(a.laplace(), 1, "Laplace transform incorrect.")

    def test_subs1(self):
        """Check subs

        """

        a = s.subs(omega)

        self.assertEqual(a.expr.is_real, True, "Lost is_real.")
        self.assertEqual2(a, omega, "Substitution fail.")

    def test_subs2(self):
        """Check subs

        """

        a1 = s.subs(omega)
        a2 = s.subs(s, omega)

        self.assertEqual(a1, a2, "Substitution fail.")

        a3 = s.subs({s: omega})
        self.assertEqual(a1, a3, "Substitution fail with dict.")

    def test_subs_const(self):
        """Check subs of a constant

        """
        a = expr('V1')
        b = expr(5 * t)
        c = a.subs({'V1': b})

        self.assertEqual(c, b, "Substitution for constant fail.")

    def test_types(self):
        """Check types

        """

        c = cexpr(10)
        self.assertEqual(type(10 + c), ConstantDomainExpression, "Not cExpr")
        self.assertEqual(type(c + 10), ConstantDomainExpression, "Not cExpr")
        self.assertEqual(type(ConstantDomainExpression(10) + c),
                         ConstantDomainExpression, "Not cExpr")
        self.assertEqual(type(c + ConstantDomainExpression(10)),
                         ConstantDomainExpression, "Not cExpr")

        self.assertEqual(type(10 + s), LaplaceDomainExpression, "Not sExpr")
        self.assertEqual(type(s + 10), LaplaceDomainExpression, "Not sExpr")
        self.assertEqual(type(ConstantDomainExpression(10) + s),
                         LaplaceDomainExpression, "Not sExpr")
        self.assertEqual(type(s + ConstantDomainExpression(10)),
                         LaplaceDomainExpression, "Not sExpr")

        self.assertEqual(type(10 + t), TimeDomainExpression, "Not tExpr")
        self.assertEqual(type(t + 10), TimeDomainExpression, "Not tExpr")
        self.assertEqual(type(ConstantDomainExpression(10) + t),
                         TimeDomainExpression, "Not tExpr")
        self.assertEqual(type(t + ConstantDomainExpression(10)),
                         TimeDomainExpression, "Not tExpr")

        v = LaplaceDomainVoltage(10)
        #self.assertEqual(type(10 + v), LaplaceDomainVoltage, "Not Vs")
        #self.assertEqual(type(v + 10), LaplaceDomainVoltage, "Not Vs")

        self.assertEqual(type(omega * t), TimeDomainExpression, "Not tExpr")
        self.assertEqual(type(t * omega), TimeDomainExpression, "Not tExpr")

        loose = state.loose_units
        state.loose_units = True
        self.assertEqual(type(LaplaceDomainExpression(10) + v),
                         LaplaceDomainVoltage, "Not Vs")
        self.assertEqual(type(v + LaplaceDomainExpression(10)),
                         LaplaceDomainVoltage, "Not Vs")
        state.loose_units = loose

    def test_evaluate(self):
        """Check evaluate

        """

        self.assertEqual(t.evaluate(10), 10.0, "Evaluate fail for scalar")
        self.assertEqual(t.evaluate((10, 20))[
                         0], 10.0, "Evaluate fail for vector")
        self.assertEqual((t * 5 + 1).evaluate((10, 20))
                         [0], 51.0, "Evaluate fail for vector")
        a = exp(t)
        self.assertEqual(a.evaluate(0j), 1, "Evaluate fail for exp(0j)")
        a = sqrt(t)
        self.assertEqual(a.evaluate(0j), 0j, "Evaluate fail for sqrt(0j)")
        self.assertEqual(a.evaluate(2j), 1 + 1j,
                         "Evaluate fail for sqrt(1+1j)")
        self.assertEqual(a.evaluate(4), 2, "Evaluate fail for sqrt(4)")

    def test_zp2k(self):
        """Test zp2k"""

        self.assertEqual(zp2tf([], [0, -1]), 1 / (s * (s + 1)), "zp2tf")

    def test_const_rms(self):
        """Test const rms"""

        c = expr(2)

        self.assertEqual(c.rms(), 2, "const rms")

    def test_assumptions(self):
        """Check assumptions

        """

        self.assertEqual(Heaviside(t).is_causal, True,
                         "Heaviside(t).is_causal")
        self.assertEqual(Heaviside(t).is_dc, False, "Heaviside(t).is_dc")
        self.assertEqual(Heaviside(t).is_ac, False, "Heaviside(t).is_ac")
        self.assertEqual(DiracDelta(t).is_causal, True,
                         "DiracDelta(t).is_causal")
        self.assertEqual(DiracDelta(t).is_dc, False, "DiracDelta(t).is_dc")
        self.assertEqual(DiracDelta(t).is_ac, False, "DiracDelta(t).is_ac")
        self.assertEqual(cos(t).is_causal, False, "cos(t).is_causal")
        self.assertEqual(cos(t).is_dc, False, "cos(t).is_dc")
        self.assertEqual(cos(t).is_ac, True, "cos(t).is_ac")

        self.assertEqual(Heaviside(t).laplace().is_causal, True,
                         "Heaviside(t).laplace().is_causal")
        self.assertEqual((Heaviside(t).laplace() +
                          DiracDelta(t).laplace()).is_causal, True)
        self.assertEqual((Heaviside(t).laplace() *
                          DiracDelta(t).laplace()).is_causal, True)
        self.assertEqual((Heaviside(t).laplace() *
                          cos(t).laplace()).is_causal, False)
        self.assertEqual(Heaviside(2 * t).is_causal, True,
                         "Heaviside(2 * t).is_causal")
        self.assertEqual(Heaviside(t - 1).is_causal, True,
                         "Heaviside(t - 1).is_causal")
        self.assertEqual(Heaviside(t + 1).is_causal, False,
                         "Heaviside(t + 1).is_causal")
        self.assertEqual(Heaviside(2 * t - 1).is_causal, True,
                         "Heaviside(2 * t - 1).is_causal")
        self.assertEqual(Heaviside(-2 * t - 1).is_causal,
                         False, "Heaviside(-2 * t - 1).is_causal")

        self.assertEqual(((Heaviside(t) - Heaviside(t - 7)) *
                         cos(t)).is_causal, True, "expr.is_causal")

    def test_has(self):
        """Test has"""

        a = Expr('3 * exp(-t) * t * a')
        self.assertEqual(a.has(3), True, "has(3)")
        self.assertEqual(a.has(4), False, "has(4)")
        self.assertEqual(a.has(t), True, "has(t)")
        self.assertEqual(a.has_symbol(t), True, "has_symbol(t)")
        self.assertEqual(a.has_symbol('a'), True, "has_symbol(a)")
        self.assertEqual(a.has_symbol('b'), False, "has_symbol(b)")

    def test_expr(self):
        """Test expr"""

        a = expr('3 * exp(-t) * t * a')
        self.assertEqual(isinstance(a, TimeDomainExpression), True, "tExpr")
        self.assertEqual(a.has(3), True, "has(3)")
        self.assertEqual(a.has(4), False, "has(4)")
        self.assertEqual(a.has(t), True, "has(t)")
        self.assertEqual(a.has_symbol(t), True, "has_symbol(t)")
        self.assertEqual(a.has_symbol('a'), True, "has_symbol(a)")
        self.assertEqual(a.has_symbol('b'), False, "has_symbol(b)")
        #self.assertEqual(a.symbols[1].is_positive, True, "a is positive")
        self.assertEqual(a.is_constant, False, "is_constant")
        self.assertEqual(a.is_unchanging, False, "is_unchanging")
        self.assertEqual(a.is_realizable, True, "is_realizable")
        # LT ignore results for t < 0
        self.assertEqual(a.is_stable, True, "is_stable")

        a = expr('3')
        self.assertEqual(a.is_constant, True, "3 is_constant")
        self.assertEqual(a.is_unchanging, True, "3 is_unchanging")

        a = expr('3 * x')
        self.assertEqual(a.is_constant, False, "3 * x is_constant")
        self.assertEqual(a.is_unchanging, True, "3 * x is_unchanging")

        self.assertEqual(isinstance(
            expr(t), TimeDomainExpression), True, "tExpr")
        self.assertEqual(isinstance(
            expr(s), LaplaceDomainExpression), True, "sExpr")
        self.assertEqual(isinstance(
            expr(f), FourierDomainExpression), True, "fExpr")
        self.assertEqual(isinstance(
            expr(omega), AngularFourierDomainExpression), True, "omegaExpr")

        self.assertEqual(isinstance(symbol('c'), Expr), True, "symbol")

        if expr('f') is not f:
            raise AssertionError('Not f')
        if expr('s') is not s:
            raise AssertionError('Not s')
        if expr('t') is not t:
            raise AssertionError('Not t')
        if expr('omega') is not omega:
            raise AssertionError('Not omega')

    def test_parallel(self):
        """Test parallel"""

        self.assertEqual(expr('4').parallel(expr('4')), 2, "parallel")

    def test_limit(self):
        """Test limit"""

        self.assertEqual(expr('4').limit(t, 0), 4, "limit")
        self.assertEqual(expr('t + 4').limit(t, 0), 4, "limit")

    def test_parameterize(self):
        """Test parameterize"""

        a = 2 / (s + 3)
        p, defs = a.parameterize()
        self.assertEqual(a, p.subs(defs), "parameterize")

        a = (s + 3) / 6
        p, defs = a.parameterize()
        self.assertEqual(a, p.subs(defs), "parameterize")

        a = (2 * s + 3) / (s + 4)
        p, defs = a.parameterize()
        self.assertEqual(a, p.subs(defs), "parameterize")

        a = (s + 3) / (s**2 + 4 * s + 4)
        p, defs = a.parameterize()
        self.assertEqual(a, p.subs(defs), "parameterize")

    def test_integrate(self):
        """Test integrate"""

        self.assertEqual(t.integrate((t, 0, t)), t**2 / 2, "integrate t")

    def test_partfrac(self):
        """Test partfrac"""

        H = expr('F(s)') / (s**2 + 3 * s + 6)

        self.assertEqual(H.partfrac(), H,  "undef partfrac")
        self.assertEqual(H.partfrac(True), H,  "undef partfrac")

        G = H * exp(-3 * s)

        self.assertEqual(G.partfrac(), G,  "undef delay partfrac")
        self.assertEqual(G.partfrac(True), G,  "undef delay partfrac")

        F = G + exp(-4 * s)

        self.assertEqual(F.partfrac(), F,  "undef delay sum partfrac")
        self.assertEqual(F.partfrac(True), F,  "undef delay sum partfrac")

        H = expr('F(s)') / (s**2 + 3 * s + 2)

        self.assertEqual(H.partfrac(), H,  "undef partfrac")
        self.assertEqual(H.partfrac(True), H,  "undef partfrac")

        G = H * exp(-3 * s)

        self.assertEqual(G.partfrac(), G,  "undef delay partfrac")
        self.assertEqual(G.partfrac(True), G,  "undef delay partfrac")

        F = G + exp(-4 * s)

        self.assertEqual(F.partfrac(), F,  "undef delay sum partfrac")
        self.assertEqual(F.partfrac(True), F,  "undef delay sum partfrac")

        H = expr('F(s)') / ((s + 2)**2)

        self.assertEqual(H.partfrac(), H,  "undef partfrac")
        self.assertEqual(H.partfrac(True), H,  "undef partfrac")

        G = H * exp(-3 * s)

        self.assertEqual(G.partfrac(), G,  "undef delay partfrac")
        self.assertEqual(G.partfrac(True), G,  "undef delay partfrac")

        F = G + exp(-4 * s)

        self.assertEqual(F.partfrac(), F,  "undef delay sum partfrac")
        self.assertEqual(F.partfrac(True), F,  "undef delay sum partfrac")

    def test_mixedfrac(self):
        """Test mixedfrac"""

        H = expr('F(s)') / (s**2 + 3 * s + 6)

        self.assertEqual(H.mixedfrac(), H,  "undef mixedfrac")

        G = H * exp(-3 * s)

        self.assertEqual(G.mixedfrac(), G,  "undef delay mixedfrac")

        F = G + exp(-4 * s)

        self.assertEqual(F.mixedfrac(), F,  "undef delay sum mixedfrac")

        H = expr('F(s)') / (s**2 + 3 * s + 2)

        self.assertEqual(H.mixedfrac(), H,  "undef mixedfrac")

        G = H * exp(-3 * s)

        self.assertEqual(G.mixedfrac(), G,  "undef delay mixedfrac")

        F = G + exp(-4 * s)

        self.assertEqual(F.mixedfrac(), F,  "undef delay sum mixedfrac")

        H = expr('F(s)') / ((s + 2)**2)

        self.assertEqual(H.mixedfrac(), H,  "undef mixedfrac")

        G = H * exp(-3 * s)

        self.assertEqual(G.mixedfrac(), G,  "undef delay mixedfrac")

        F = G + exp(-4 * s)

        self.assertEqual(F.mixedfrac(), F,  "undef delay sum mixedfrac")

    def force_causal(self):
        """Test force_causal"""

        X = 1 / (1 + s)
        x = X(t)

        self.assertEqual(x.force_causal(), exp(-t) *
                         Heaviside(t), "force causal")

        x = X(t, causal=True)
        self.assertEqual(x.force_causal(), exp(-t) *
                         Heaviside(t), "force causal if already causal")

    def strip_condition(self):
        """Test strip_condition"""

        X = 1 / (1 + s)
        x = X(t)

        self.assertEqual(x.strip_condition(), exp(-t), "remove condition")

        x = X(t, condition=True)
        self.assertEqual(x.strip_condition(), exp(-t),
                         "remove condition if causal")

    def test_val(self):
        """Test val"""

        a = expr('3 / 4')

        self.assertEqual(a.val, 0.75, "val")
        self.assertEqual(a.evalf(), 0.75, "evalf")

    def test_cval(self):
        """Test cval"""

        a = 3j
        e = expr(a)

        self.assertEqual(e.cval, a, 'cval')

    def test_fval(self):
        """Test fval"""

        a = 3
        e = expr(a)

        self.assertEqual(e.fval, a, 'fval')

    def test_equality(self):
        """Test equality"""

        e = Eq(expr('x(t)'), cos(3 * t))
        self.assertEqual(e.is_ac, True, 'is_ac')
        self.assertEqual(e.is_dc, False, 'is_dc')
        self.assertEqual(e.is_causal, False, 'is_causal')

    def test_phase(self):
        """Test phase"""

        c = expr(2).phase
        e = exp(c)

        self.assertEqual(c.part, 'phase', "is_phase")
        self.assertEqual(e.part, '', "exp is_phase")

    def test_parts(self):
        """Test parts"""

        self.assertEqual(expr(1).real.is_real_part, True, "is_real_part")
        self.assertEqual(expr(1).imag.is_imag_part, True, "is_imag_part")
        self.assertEqual(expr(1).magnitude.is_magnitude, True, "is_magnitude")
        self.assertEqual(expr(1).phase.is_phase, True, "is_phase")
        self.assertEqual(expr(1).phase_degrees.is_phase_degrees,
                         True, "is_phase_degrees")
        self.assertEqual(expr(1).dB.is_dB, True, "is_dB")

    def test_dB(self):
        """Test dB"""

        v = voltage(10)
        self.assertEqual(v.dB, 20, "voltage(10).dB")
        self.assertEqual((v**2).dB, 20, "(voltage(10)**2).dB")

    def test_poles(self):
        """Test poles"""

        # TODO, need to handle different orderings of complex conjugate pairs

        # H = 1 / ((s + 5j)**3 * (s-5j)**2 * (s+2))
        # pairs, singles = H.poles(pairs=True)
        # epairs = D({(5 * j, -5 * j): 2})
        # esingles = D({-2: 1, -5 * j: 1})
        # self.assertEqual(set(pairs.keys()), set(epairs.keys()), "pole pairs")
        # self.assertEqual(set(singles.keys()), set(esingles.keys()), "pole singles")

        # H = 1 / ((s + 5j)**2 * (s-5j)**3 * (s+2))
        # pairs, singles = H.poles(pairs=True)
        # epairs = D({(5 * j, -5 * j): 2})
        # esingles = D({-2: 1, 5 * j: 1})
        # self.assertEqual(set(pairs.keys()), set(epairs.keys()), "pole pairs")
        # self.assertEqual(set(singles.keys()), set(esingles.keys()), "pole singles")

        # H = 1 / ((s + 5j)**2 * (s-5j)**2 * (s+2))
        # pairs, singles = H.poles(pairs=True)
        # epairs = D({(5 * j, -5 * j): 2})
        # esingles = D({-2: 1})
        # self.assertEqual(set(pairs.keys()), set(epairs.keys()), "pole pairs")
        # self.assertEqual(set(singles.keys()), set(esingles.keys()), "pole singles")

    def test_auto_quantity(self):

        a = 4 * t * volts
        self.assertEqual(a.as_voltage().is_voltage, True, 'is_voltage')

        a = 4 * t * amperes
        self.assertEqual(a.as_current().is_current, True, 'is_current')

    def test_pole_zero_units(self):
        """Test poles and zeros units"""

        H = s / (s + 3)
        self.assertEqual(str(list(H.poles())[0].units), 'rad/s', "pole units")
        self.assertEqual(str(list(H.zeros())[0].units), 'rad/s', "zero units")

    def test_uExpr1(self):
        """Check uExpr

        """
        a = uexpr('x * 3 + 4', 'x')
        b = expr('x * 3 + 4', 'x')

        x = symbol('x').sympy

        c = b * b

        self.assertEqual(a, b, 'uexpr')
        self.assertEqual(a.var, x, 'uexpr.var')
        self.assertEqual(c.var, x, '(uexpr * uexpr).var')
