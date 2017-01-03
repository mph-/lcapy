from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def assertEqual2(self, ans1, ans2, comment):

        try:
            self.assertEqual(ans1, ans2, comment)
        except AssertionError as e:
            pprint(ans1)
            pprint(ans2)
            raise AssertionError(e)

    def test_rmul(self):

        self.assertEqual2((j * omega).__class__,  omegaExpr, "j * omega fail.")

    def test_sExpr1(self):
        """Lcapy: check sExpr1

        """
        a = sExpr('(s+2)/(s-2)')
        self.assertEqual2(a.N, sExpr('s+2'), "N incorrect.")
        self.assertEqual2(a.D, sExpr('s-2'), "D incorrect.")

        self.assertEqual2(sorted(a.poles()), [2], "poles incorrect.")
        self.assertEqual2(sorted(a.zeros()), [-2], "zeros incorrect.")

        self.assertEqual2(a.partfrac(), 1 + 4 / (s - 2), "partfrac incorrect.")

        self.assertEqual(a.evaluate(1), -3.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate(-2), 0.0, "scalar evaluate incorrect.")

        self.assertEqual(a.inverse_laplace(causal=True), 4 * exp(2 * t) * H(t) + DiracDelta(t), "inverse Laplace incorrect.")

        aw = a(j * omega)
        self.assertEqual2(aw.real, (omega**2 - 4) / (omega**2 + 4), "real part incorrect.")
        self.assertEqual2(aw.imag, -4 * omega / (omega**2 + 4), "imag part incorrect.")


    def test_sExpr2(self):
        """Lcapy: check sExpr2

        """
        a = (s + 2) * (s + 3) / (s - 2)
        self.assertEqual2(a.N, (s + 2) * (s + 3), "N incorrect.")
        self.assertEqual2(a.D, s - 2, "D incorrect.")

        self.assertEqual2(sorted(a.poles()), [2], "poles incorrect.")
        self.assertEqual2(sorted(a.zeros()), sorted([-3.0, -2.0]), "zeros incorrect.")

        self.assertEqual2(
            a.partfrac(), s + 7 + 20 / (s - 2), "partfrac incorrect.")
        self.assertEqual2(
            a.mixedfrac(), s + 7 + 20 / (s - 2), "mixedfrac incorrect.")
        self.assertEqual2(
            a.general(), (s**2 + 5 * s + 6) / (s - 2), "general incorrect.")
        self.assertEqual2(
            a.canonical(), (s**2 + 5 * s + 6) / (s - 2), "canonical incorrect.")

        self.assertEqual(a.inverse_laplace(causal=True), 20 * exp(2 * t) * H(t) + 7 * DiracDelta(t) + DiracDelta(t, 1), "inverse Laplace incorrect.")

    def test_sExpr3(self):
        """Lcapy: check sExpr3

        """
        a = (s**2 + 5 * s + 6) / (s - 2)
        self.assertEqual2(a.N, s ** 2 + 5 * s + 6, "N incorrect.")
        self.assertEqual2(a.D, s - 2, "D incorrect.")

        self.assertEqual2(sorted(a.poles()), [2], "poles incorrect.")
        self.assertEqual2(sorted(a.zeros()), sorted([-3, -2]), "zeros incorrect.")

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

        self.assertEqual(a.inverse_laplace(causal=True), 20 * exp(2 * t) * H(t) + 7 * DiracDelta(t) + DiracDelta(t, 1), "inverse Laplace incorrect.")

    def test_sExpr4(self):
        """Lcapy: check sExpr4

        """
        a = 1 / ((s - j) * (s + j))
        self.assertEqual2(a.N, 1, "N incorrect.")
        self.assertEqual2(a.D, (s - j) * (s + j), "D incorrect.")

        self.assertEqual2(set(a.poles()), set([-j, j]), "poles incorrect.")
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

        self.assertEqual(a.inverse_laplace(causal=True), sin(t) * H(t), "inverse Laplace incorrect.")

    def test_sExpr5(self):
        """Lcapy: check sExpr5

        """
        a = 1 / ((s - 1j) * (s + 1j))
        self.assertEqual2(a.N, 1, "N incorrect.")
        self.assertEqual2(a.D, (s - j) * (s + j), "D incorrect.")

        self.assertEqual2(set(a.poles()), set([-j, j]), "poles incorrect.")
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

        self.assertEqual(a.inverse_laplace(causal=True), sin(t) * H(t), "inverse Laplace incorrect.")

    def test_wExpr1(self):
        """Lcapy: check wExpr1

        """

        A = (j * omega + 3) / (j * omega - 4)
        self.assertEqual2(A.N, j * omega + 3, "N incorrect.")
        self.assertEqual2(A.D, j * omega - 4, "D incorrect.")
        self.assertEqual2(A.real, (omega**2 - 12) / (omega**2 + 16), "real incorrect.")
        self.assertEqual2(A.imag, -7 * omega / (omega**2 + 16), "imag incorrect.")

    def test_tExpr1(self):
        """Lcapy: check tExpr1

        """
        a = t**2

        self.assertEqual(a.evaluate(2), 4.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate((2, 3))[1], 9.0, "vector evaluate incorrect.")

    def test_step(self):
        """Lcapy: check step

        """
        a = u(t)

        self.assertEqual(a.evaluate(2), 1.0, "scalar evaluate incorrect.")
        self.assertEqual(a.evaluate((-2, 2))[1], 1.0, "vector evaluate incorrect.")
        self.assertEqual(a.evaluate((-2, 2))[0], 0.0, "vector evaluate incorrect.")
        self.assertEqual(a.laplace(), 1 / s, "Laplace transform incorrect.")

    def test_delta(self):
        """Lcapy: check delta

        """
        a = delta(t)

        self.assertEqual(a.evaluate(2), 0, "scalar evaluate incorrect.")
        self.assertEqual(a.laplace(), 1, "Laplace transform incorrect.")

    def test_jomega(self):
        """Lcapy: check jomega

        """

        a = sExpr('s+2')
        b = a(j * omega)

        self.assertEqual2(b, j * omega + 2, "Substitution failed.")
        self.assertEqual2(a.jomega, j * omega + 2, "jomega failed.")

    def test_subs1(self):
        """Lcapy: check subs

        """

        a = s(omega)

        self.assertEqual(a.expr.is_real, True, "Lost is_real.")
        self.assertEqual2(a, omega, "Substitution fail.")

    def test_subs2(self):
        """Lcapy: check subs

        """

        a1 = s(omega)
        a2 = s.subs(s, omega)

        self.assertEqual(a1, a2, "Substitution fail.")

        a3 = s.subs({s: omega})
        self.assertEqual(a1, a3, "Substitution fail with dict.")


    def test_types(self):
        """Lcapy: check types

        """

        c = cExpr(10)
        self.assertEqual(type(10 + c), cExpr, "Not cExpr")
        self.assertEqual(type(c + 10), cExpr, "Not cExpr")
        self.assertEqual(type(cExpr(10) + c), cExpr, "Not cExpr")
        self.assertEqual(type(c + cExpr(10)), cExpr, "Not cExpr")

        self.assertEqual(type(10 + s), sExpr, "Not sExpr")
        self.assertEqual(type(s + 10), sExpr, "Not sExpr")
        self.assertEqual(type(cExpr(10) + s), sExpr, "Not sExpr")
        self.assertEqual(type(s + cExpr(10)), sExpr, "Not sExpr")

        self.assertEqual(type(10 + t), tExpr, "Not tExpr")
        self.assertEqual(type(t + 10), tExpr, "Not tExpr")
        self.assertEqual(type(cExpr(10) + t), tExpr, "Not tExpr")
        self.assertEqual(type(t + cExpr(10)), tExpr, "Not tExpr")

        v = Vs(10)
        self.assertEqual(type(10 + v), Vs, "Not Vs")
        self.assertEqual(type(v + 10), Vs, "Not Vs")
        self.assertEqual(type(sExpr(10) + v), Vs, "Not Vs")
        self.assertEqual(type(v + sExpr(10)), Vs, "Not Vs")

        self.assertEqual(type(omega * t), tExpr, "Not tExpr")
        self.assertEqual(type(t * omega), tExpr, "Not tExpr")                

    def test_evaluate(self):
        """Lcapy: check evaluate

        """

        self.assertEqual(t.evaluate(10), 10.0, "Evaluate fail for scalar")
        self.assertEqual(t.evaluate((10, 20))[0], 10.0, "Evaluate fail for vector")
        self.assertEqual((t * 5 + 1).evaluate((10, 20))[0], 51.0, "Evaluate fail for vector")

    def test_laplace(self):

        self.assertEqual(Heaviside(t).laplace(), 1 / s, "Heaviside(t)")
        self.assertEqual(DiracDelta(t).laplace(), 1, "DiracDelta(t)")

    def test_zp2k(self):

        self.assertEqual(zp2tf([], [0, -1]), 1 / (s * (s + 1)), "zp2tf")

    def test_fourier(self):

        self.assertEqual((t * 0 + 1).fourier(), DiracDelta(f))
        self.assertEqual((t * 0 + 1).fourier().inverse_fourier(), 1)
        self.assertEqual(t.fourier(), 2 * j * pi * DiracDelta(f, 1))
        self.assertEqual(2 * cos(2 * pi * t).fourier(),
                         DiracDelta(f - 1) + DiracDelta(f + 1))
        self.assertEqual(2 * sin(2 * pi * t).fourier(),
                         -j * DiracDelta(f - 1) + j * DiracDelta(f + 1))
        self.assertEqual(exp(j * 2 * pi * t).fourier(), DiracDelta(f - 1))
        self.assertEqual(exp(j * 2 * pi * t).fourier().inverse_fourier(),
                         exp(j * 2 * pi * t))
        self.assertEqual(exp(-j * 2 * pi * t).fourier(), DiracDelta(f + 1))
        self.assertEqual(exp(-j * 2 * pi * t).fourier().inverse_fourier(),
                         exp(-j * 2 * pi * t))

    def test_inverse_fourier(self):

        self.assertEqual((1 / (s + 1))(j * omega).inverse_fourier(), exp(-t) * Heaviside(t))
        self.assertEqual((1 / (s + 1))(j * omega)(2 * pi * f).inverse_fourier(), exp(-t) * Heaviside(t))

    def test_rms(self):

        self.assertEqual(Vconst(2).rms(), Vt(2))
        self.assertEqual(Vphasor(2).rms(), Vt(1))

    def test_assumptions(self):
        """Lcapy: check assumptions

        """

        self.assertEqual(Heaviside(t).is_causal, True, "Heaviside(t).is_causal")
        self.assertEqual(Heaviside(t).is_dc, False, "Heaviside(t).is_dc")
        self.assertEqual(Heaviside(t).is_ac, False, "Heaviside(t).is_ac")
        self.assertEqual(DiracDelta(t).is_causal, True, "DiracDelta(t).is_causal")
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
                          cos(t).laplace()).is_causal, True)
