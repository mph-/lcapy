from lcapy import *
from lcapy.discretetime import *
import unittest
import sympy as sym


class LcapyZTransforms(unittest.TestCase):

    """Unit tests for lcapy ztransforms

    """

    def test_ztransform(self):

        self.assertEqual(nexpr(1).ZT(), z / (z - 1), "1")
        self.assertEqual(nexpr(3).ZT(), 3 * z / (z - 1), "3")
        self.assertEqual(Heaviside(n).ZT(), z / (z - 1), "Heaviside(n)")
        self.assertEqual(3 * Heaviside(n).ZT(), 3 *
                         z / (z - 1), "3 * Heaviside(n)")
        self.assertEqual(Heaviside(n + 1).ZT(), z **
                         2 / (z - 1), "Heaviside(n + 1)")
        self.assertEqual(Heaviside(n - 1).ZT(), 1 /
                         (z - 1), "Heaviside(n - 1)")
        self.assertEqual(unitimpulse(n).ZT(), 1, "unitimpulse(n)")
        self.assertEqual(3 * unitimpulse(n).ZT(), 3, "unitimpulse(n)")
        self.assertEqual(unitimpulse(n - 1).ZT(), 1 / z, "unitimpulse(n - 1)")
        self.assertEqual(unitimpulse(n - 2).ZT(), 1 /
                         z**2, "unitimpulse(n - 2)")
        self.assertEqual(3 * unitimpulse(n - 2).ZT(), 3 /
                         z**2, "3 * unitimpulse(n - 2)")
        self.assertEqual(3 * unitimpulse(n + 2).ZT(), 3 *
                         z**2, "3 * unitimpulse(n + 2)")
        self.assertEqual(expr('v(n)').ZT(), expr('V(z)'), "v(n)")
        self.assertEqual(expr('v(n / 3)').ZT(), expr('V(z**3)'), "v(n/3)")
        self.assertEqual(expr('3 * v(n / 3)').ZT(),
                         expr('3 * V(z**3)'), "3 * v(n/3)")
        self.assertEqual(expr('v(n-3)').ZT(), expr('V(z) / z**3'), "v(n - 3)")
        self.assertEqual(nexpr('x(n)').ZT(), zexpr('X(z)'), "x(n)")

        self.assertEqual((0.1**n).ZT(), 10 * z / (10 * z - 1), "0.1**n")
        self.assertEqual(((0.1**n) * u(n)).ZT(), 10 *
                         z / (10 * z - 1), "0.1**n")

        self.assertEqual(exp(3 * n).ZT(), z / (z - exp(3)), "exp(3 * n)")
        self.assertEqual((exp(3 * n) * exp(2 * j * n)).ZT(), z /
                         (z - exp(3 + 2 * j)), "exp(3 * n) * exp(2 * j * n)")
        self.assertEqual(cos(3 * n).ZT(), (z**2 - z * cos(3)) /
                         (z**2 - 2 * z * cos(3) + 1), "cos(3 * n)")
        self.assertEqual(sin(3 * n).ZT(), (z * sin(3)) /
                         (z**2 - 2 * z * cos(3) + 1), "sin(3 * n)")

        e = 5**n
        self.assertEqual(e(z)(n).remove_condition(), e, "5**n")
        e = 5**-n
        self.assertEqual(e(z)(n).remove_condition(), e, "5**-n")

        e = n * 5**n
        self.assertEqual(e(z)(n).remove_condition(), e, "n * 5**n")
        e = n * 5**-n
        self.assertEqual(e(z)(n).remove_condition(), e, "n * 5**-n")

        e = n * exp(3 * n)
        self.assertEqual(e(z)(n).remove_condition(), e, "n * exp(3 * n)")

        # Added by Juergen Weizenecker
        self.assertEqual((n).ZT(), z/(z - 1)**2, "n")
        self.assertEqual((n**2).ZT(), z*(z + 1)/(z - 1)**3, "n**2")
        self.assertEqual((8*n**3).ZT(), 8*z*(3*z*(z + 1) -
                         (z - 1)*(2*z + 1))/(z - 1)**4, "8*n**3")
        self.assertEqual(((2*n + 3)**2).ZT(), z *
                         (16*z + 9*(z - 1)**2 - 8)/(z - 1)**3, "(2*n + 3)**2")
        self.assertEqual((2**(3*n + 4)).ZT(), 16*z/(z - 8), "2**(3*n + 4)")
        self.assertEqual((2**(3*n + 4)*n).ZT(), 128 *
                         z/(z - 8)**2, "2**(3*n + 4)*n")
        self.assertEqual((3**(4*n + 5)*n**2).ZT(), 19683*z *
                         (z + 81)/(z - 81)**3, "3**(4*n + 5)*n**2")
        self.assertEqual((sin(2*n + 3)).ZT(), z*(z*sin(3) -
                         sin(1))/(z**2 - 2*z*cos(2) + 1), "sin(2*n + 3)")
        self.assertEqual((cos(2*n + 3)).ZT(), z*(z*cos(3) -
                         cos(1))/(z**2 - 2*z*cos(2) + 1), "cos(2*n + 3)")
        self.assertEqual((n*sin(2*n + 3)).ZT(), z*(z**2*sin(5) - 2*z *
                         sin(3) + sin(1))/(z**2 - 2*z*cos(2) + 1)**2, "n*sin(2*n + 3)")
        self.assertEqual((n*cos(2*n + 3)).ZT(), z*(z**2*cos(5) - 2*z *
                         cos(3) + cos(1))/(z**2 - 2*z*cos(2) + 1)**2, "n*cos(2*n + 3)")
        self.assertEqual((2**(3*n + 4)*sin(5*n + 6)).ZT(), 16*z*(z*sin(6) -
                         8*sin(1))/(z**2 - 16*z*cos(5) + 64), "2**(3*n + 4)*sin(5*n + 6)")
        self.assertEqual((2**(3*n + 4)*cos(5*n + 6)).ZT(), 16*z*(z*cos(6) -
                         8*cos(1))/(z**2 - 16*z*cos(5) + 64), "2**(3*n + 4)*cos(5*n + 6)")
        self.assertEqual((2**(3*n + 4)*n*sin(5*n + 6)).ZT(), 128*z*(z**2*sin(11) - 16*z*sin(
            6) + 64*sin(1))/(z**2 - 16*z*cos(5) + 64)**2, "2**(3*n + 4)*n*sin(5*n + 6)")
        self.assertEqual((2**(3*n + 4)*n*cos(5*n + 6)).ZT(), 128*z*(z**2*cos(11) - 16*z*cos(
            6) + 64*cos(1))/(z**2 - 16*z*cos(5) + 64)**2, "2**(3*n + 4)*n*cos(5*n + 6)")
        self.assertEqual((3**(4*n + 5)*n**2*sin(6*n + 7)).ZT(), 19683*z*(4*z*(z - 81*cos(6))*(z**2*sin(13) - 162*z*sin(7) + 6561*sin(1)) - 3*(
            z**2 - 162*z*cos(6) + 6561)*(z**2*sin(13) - 108*z*sin(7) + 2187*sin(1)))/(z**2 - 162*z*cos(6) + 6561)**3, "3**(4*n + 5)*n**2*sin(6*n + 7)")
        self.assertEqual((3**(4*n + 5)*n**2*cos(6*n + 7)).ZT(), 19683*z*(4*z*(z - 81*cos(6))*(z**2*cos(13) - 162*z*cos(7) + 6561*cos(1)) - 3*(
            z**2 - 162*z*cos(6) + 6561)*(z**2*cos(13) - 108*z*cos(7) + 2187*cos(1)))/(z**2 - 162*z*cos(6) + 6561)**3, "3**(4*n + 5)*n**2*cos(6*n + 7)")
        self.assertEqual((n*exp(3*n + 4)).ZT(), z*exp(7) /
                         (z - exp(3))**2, "n*exp(3*n + 4)")
        self.assertEqual((n**2*exp(3*n + 4)).ZT(), z*(z + exp(3))
                         * exp(7)/(z - exp(3))**3, "n**2*exp(3*n + 4)")
        self.assertEqual((exp(2*n + 3)*sin(5*n + 4)).ZT(), z*(z*sin(4) + exp(2)*sin(1))
                         * exp(3)/(z**2 - 2*z*exp(2)*cos(5) + exp(4)), "exp(2*n + 3)*sin(5*n + 4)")
        self.assertEqual((exp(3*n + 2)*cos(5*n + 4)).ZT(), z*(z*cos(4) - exp(3)*cos(1))
                         * exp(2)/(z**2 - 2*z*exp(3)*cos(5) + exp(6)), "exp(3*n + 2)*cos(5*n + 4)")
        self.assertEqual((n*exp(2*n + 3)*sin(5*n + 4)).ZT(), -z*(-z**2*sin(9) + 2*z*exp(2)*sin(4) + exp(
            4)*sin(1))*exp(5)/(z**2 - 2*z*exp(2)*cos(5) + exp(4))**2, "n*exp(2*n + 3)*sin(5*n + 4)")
        self.assertEqual((n*exp(3*n + 2)*cos(5*n + 4)).ZT(), z*(z**2*cos(9) - 2*z*exp(3)*cos(4) + exp(
            6)*cos(1))*exp(5)/(z**2 - 2*z*exp(3)*cos(5) + exp(6))**2, "n*exp(3*n + 2)*cos(5*n + 4)")
        self.assertEqual((n**2*exp(2*n + 3)*sin(5*n + 4)).ZT(), z*(z**4*sin(9) + z**3*exp(2)*sin(14) - 3*z**3*exp(2)*sin(4) - 3*z**2*exp(4)*sin(1) - 3*z**2*exp(
            4)*sin(9) + 3*z*exp(6)*sin(4) + z*exp(6)*sin(6) + exp(8)*sin(1))*exp(5)/(z**2 - 2*z*exp(2)*cos(5) + exp(4))**3, "n**2*exp(2*n + 3)*sin(5*n + 4)")
        self.assertEqual((exp(-2*n - 5)*exp(4*n + 2)).ZT(),
                         z*exp(-3)/(z - exp(2)), "exp(-2*n - 5)*exp(4*n + 2)")
        self.assertEqual(((2 + 3*j)**n).ZT(), -z /
                         (-z + 2 + 3*j), "(2 + 3*j)**n")
        self.assertEqual((exp(n*(3 + 2*j))).ZT(), z /
                         (z - exp(3 + 2*j)), "exp(n*(3 + 2*j))")
        self.assertEqual((n*(2 + 3*j)**n).ZT(), z*(2 + 3*j) /
                         (-z + 2 + 3*j)**2, "n*(2 + 3*j)**n")
        self.assertEqual((n*exp(n*(3 + 2*j))).ZT(), z*exp(3 + 2*j) /
                         (z - exp(3 + 2*j))**2, "n*exp(n*(3 + 2*j))")
        self.assertEqual((n**2*(2 + 3*j)**n).ZT(), -z*(2 + 3*j)
                         * (z + 2 + 3*j)/(-z + 2 + 3*j)**3, "n**2*(2 + 3*j)**n")
        self.assertEqual((n**2*exp(n*(3 + 2*j))).ZT(), z*(z + exp(3 + 2*j))
                         * exp(3 + 2*j)/(z - exp(3 + 2*j))**3, "n**2*exp(n*(3 + 2*j))")
        self.assertEqual((2**n*exp(3*j*n)).ZT(), z /
                         (z - 2*exp(3*j)), "2**n*exp(3*j*n)")
        self.assertEqual((2**n*n*exp(3*j*n)).ZT(), 2*z*exp(3*j) /
                         (z - 2*exp(3*j))**2, "2**n*n*exp(3*j*n)")
        self.assertEqual((3**n*n**2*exp(4*j*n)).ZT(), 3*z*(z + 3*exp(4*j))
                         * exp(4*j)/(z - 3*exp(4*j))**3, "3**n*n**2*exp(4*j*n)")

    def test_inverse_ztransform(self):

        self.assertEqual(zexpr(1).IZT(causal=True), unitimpulse(n), "1")
        self.assertEqual((z**-1).IZT(causal=True), unitimpulse(n-1), "z**-1")
        self.assertEqual((z**-2).IZT(causal=True), unitimpulse(n-2), "z**-2")
        self.assertEqual(zexpr('1 / (1 - a * z ** -1)').IZT(causal=True),
                         nexpr('a**n * u(n)'), "1 / (1 - a * z)")
        self.assertEqual(zexpr('X(z)').IZT(causal=True), nexpr('x(n)'), "X(z)")

        # Added by Juergen Weizenecker
        self.assertEqual((z*(z**2 + 2*z + 1/4)/(2*(z - 1/2)**4)).IZT(causal=True,
                         pairs=False), 2**(-n)*n**3*u(n), "2**(-n)*n**3*u(n)")
        self.assertEqual((-z*(-2*z + 2 + exp(2))/((z - 2)*(z - exp(2)))).IZT(causal=True, pairs=False), -
                         (-2**n - exp(2*n))*unitstep(n), "-z*(-2*z + 2 + exp(2))/((z - 2)*(z - exp(2)))")
        self.assertEqual((z*(1 + 2*j)*(3*z*(z + 1 + 2*j) + (-z + 1 + 2*j)*(2*z + 1 + 2*j))/(-z + 1 + 2*j)**4).IZT(causal=True, pairs=False), (1 + 2*j)*(1 + 2*j)**n*(n*(-18 + 24*j)*(n - 2)*(
            n - 1)/(6*(1 + 2*j)**3) + n*(6 + 12*j)*(n - 1)/(2*(1 + 2*j)**2) + n*(1 - 2*j)/5)*unitstep(n), "z*(1 + 2*j)*(3*z*(z + 1 + 2*j) + (-z + 1 + 2*j)*(2*z + 1 + 2*j))/(-z + 1 + 2*j)**4")
        self.assertEqual((z*(6*z + 5 - exp(3 + 4*j))/((z + 1)*(z - exp(3 + 4*j)))).IZT(causal=True, pairs=False),
                         ((-1)**n + 5*exp(3 + 4*j)**n)*unitstep(n), "z*(6*z + 5 - exp(3 + 4*j))/((z + 1)*(z - exp(3 + 4*j)))")
        self.assertEqual((z*((z - 6)**4*(z + 1) + (z - 1)**3*(18*z*(z + 6) - 12*(z - 6)*(z + 3)))/((z - 6)**4*(z - 1)**3)).IZT(causal=True, pairs=False), (6**n*(
            n*(n - 2)*(n - 1) + 3*n*(n - 1) + n) + n*(n - 1) + n)*unitstep(n), "z*((z - 6)**4*(z + 1) + (z - 1)**3*(18*z*(z + 6) - 12*(z - 6)*(z + 3)))/((z - 6)**4*(z - 1)**3)")
        self.assertEqual((4*z*(18*z - 5)/((2*z - 1)*(4*z - 1))).IZT(causal=True, pairs=False),
                         4*(4**(-n)/4 + 2*2**(-n))*unitstep(n), "4*z*(18*z - 5)/((2*z - 1)*(4*z - 1))")
        self.assertEqual((z*(z - 3*cos(2))/(z**2 - 6*z*cos(2) + 9)).IZT(causal=True, pairs=True),
                         3**n*cos(2*n)*unitstep(n), "z*(z - 3*cos(2))/(z**2 - 6*z*cos(2) + 9)")
        self.assertEqual((-3*z*sin(2)/(z**2 + 6*z*cos(2) + 9)).IZT(causal=True, pairs=True), -
                         3**n*sin(n*(-2 + pi))*unitstep(n), "-3*z*sin(2)/(z**2 + 6*z*cos(2) + 9)")
        # This fails with ubuntu-18.04 python < 3.9
        # self.assertEqual((3*z*(z**4*cos(2) - 6*z**3 - 6*z**3*sin(2)**2 + 54*z*sin(2)**2 + 54*z - 81*cos(2))/(z**2 - 6*z*cos(2) + 9)**3).IZT(causal=True, pairs=True), 3**n*n**2*cos(2*n)*unitstep(n), "3*z*(z**4*cos(2) - 6*z**3 - 6*z**3*sin(2)**2 + 54*z*sin(2)**2 + 54*z - 81*cos(2))/(z**2 - 6*z*cos(2) + 9)**3")

        self.assertEqual((z*(z**2 - exp(4))*exp(5)*sin(2)/(z**2 - 2*z*exp(2)*cos(2) + exp(4))**2).IZT(causal=True, pairs=True),
                         n*exp(3)*exp(2*n)*sin(2*n)*unitstep(n),
                         "z*(z**2 - exp(4))*exp(5)*sin(2)/(z**2 - 2*z*exp(2)*cos(2) + exp(4))**2")
        self.assertEqual((3*z*(z**2*sin(6) - 6*z*sin(4) + 9*sin(2))/(z**2 - 6*z*cos(2) + 9)**2).IZT(causal=True, pairs=True), 3*3**n*(
            n*sin(2*n)*cos(4)/3 + n*sin(4)*cos(2*n)/3)*unitstep(n), "3*z*(z**2*sin(6) - 6*z*sin(4) + 9*sin(2))/(z**2 - 6*z*cos(2) + 9)**2")
        self.assertEqual((2*z*exp(3*j)/(z - 2*exp(3*j))**2).IZT(causal=True, pairs=False),
                         n*(2*exp(3*j))**n*unitstep(n), "2*z*exp(3*j)/(z - 2*exp(3*j))**2")

    def test_misc(self):

        x = 3 * unitimpulse(n - 2)
        y1 = x.differentiate()
        y2 = x.ZT().ndifferentiate().IZT()
        self.assertEqual(y1, y2, "differentiate")

        x = 3 * unitimpulse(n - 2)
        y1 = x.integrate()
        y2 = x.ZT().nintegrate().IZT()
        # Need to teach that these are equal
        #self.assertEqual(y1, y2, "integrate")

    def test_convolution(self):

        a = expr('a(n)')
        b = expr('b(n)')
        c = (a.ZT() * b.ZT()).IZT(causal=True)
        d = expr('Sum(a(-m + n)*b(m), (m, 0, n))')

        self.assertEqual(c, d, "convolution")
