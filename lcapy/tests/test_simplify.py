from lcapy import *
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

    def test_simplify_sin_cos(self):
        """Test simplify_sin_cos"""

        a = expr('3 * cos(theta) + 4 * sin(theta)')
        b = expr('5 * cos(theta - atan(4 / 3))')

        # Note, simplify converts the results to the same form.
        self.assertEqual(a.simplify_sin_cos(
            as_cos=True).simplify(), b, "cos form")
        self.assertEqual(a.simplify_sin_cos(
            as_sin=True).simplify(), b, "sin form")
        self.assertEqual(a.simplify_sin_cos().simplify(), b, "default form")

        a = 3 * sin(2 * t) + 4 * cos(2 * t) + 5 * sin(4 * t) + \
            12 * cos(4 * t) + cos(2 * t + 1) + 3
        b = 5 * cos(2 * t - atan('3 / 4')) + 13 * \
            cos(4 * t - atan('5 / 12')) + cos(2 * t + 1) + 3

        self.assertEqual(a.simplify_sin_cos().simplify(), b, "default")

    def test_simplify_heaviside(self):
        """Test simplify_heaviside"""

        a = Heaviside(4 * t - 2)
        b = a.simplify_heaviside()
        self.assertEqual(b, Heaviside(t - 0.5), "simplify_heaviside scale")

        a = Heaviside(t) * Heaviside(t)
        self.assertEqual(a.simplify_heaviside(), Heaviside(t),
                         "simplify_heaviside product")

    def test_simplify_rect(self):
        """Test simplify_rect"""

        a = rect(t) * rect(t)
        self.assertEqual(a.simplify_rect(), rect(t), "simplify_rect product")

    def test_simplify_dirac_delta(self):
        """Test simplify_dirac_delta"""

        a = DiracDelta(t) * 'g(t)'
        self.assertEqual(a.simplify_dirac_delta(), DiracDelta(
            t) * 'g(0)', "simplify_dirac_delta product")

        a = DiracDelta(3 * t)
        self.assertEqual(a.simplify_dirac_delta(), DiracDelta(
            t) / 3, "simplify_dirac_delta scale")

    def test_simplify_conjugates(self):
        """Test simplify_conjugates"""

        a = exp(j * 2) + exp(-j * 2) + exp(j * 3)
        self.assertEqual(a.simplify_conjugates(), 2 * cos(2) + exp(j * 3), "")

    def test_expand_hyperbolic_trig(self):
        """Test expand_hyperbolic_trig"""

        a = cosh(2) + sinh(2)
        self.assertEqual(a.expand_hyperbolic_trig(), 2 * exp(2), "")

        a = cosh(2) - sinh(2)
        self.assertEqual(a.expand_hyperbolic_trig(), 2 * exp(-2), "")
