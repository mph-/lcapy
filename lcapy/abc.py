"""This module exports all latin and greek letters as Symbols, so you can
conveniently do

    >>> from lcapy.abc import x, y

instead of the slightly more clunky-looking

    >>> from lcapy import symbols
    >>> x, y = symbols('x y')

Unlike sympy.abc, this does not define symbols that clash with Lcapy's
domain variables.

"""

from lcapy.expr import symbols

a, b, c, d, e, g, h, i = symbols('a, b, c, d, e, g, h, i')
l, m, o, p, q, r = symbols('l, m, o, p, q, r')
u, v, w, x, y = symbols('u, v, w, x, y')

A, B, C, D, E, G, H, I, J = symbols('A, B, C, D, E, G, H, I, J')
K, L, M, N, O, P, Q, R, S, T = symbols('K, L, M, N, O, P, Q, R, S, T')
U, V, W, X, Y, Z = symbols('U, V, W, X, Y, Z')

alpha, beta, gamma, delta = symbols('alpha, beta, gamma, delta')
epsilon, zeta, eta, theta = symbols('epsilon, zeta, eta, theta')
iota, kappa, lamda, mu = symbols('iota, kappa, lamda, mu')
nu, xi, omicron, pi = symbols('nu, xi, omicron, pi')
rho, sigma, tau, upsilon = symbols('rho, sigma, tau, upsilon')
phi, chi, psi = symbols('phi, chi, psi')
