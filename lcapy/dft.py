"""This module provides support for discrete Fourier transforms.

It calculates the discrete Fourier transform using:

   X(k) = \sum_{n=0}^{N-1} x(n) e^{-j * 2 * \pi * n * k / N}

Copyright 2020--2021 Michael Hayes, UCECE
Juergen Weizenecker, HKA
"""

import sympy as sym
from .transformer import BilateralForwardTransformer
from .sym import sympify, AppliedUndef, j, pi, miscsymbol
from .extrafunctions import UnitImpulse, UnitStep
from .utils import factor_const, scale_shift
from .matrix import Matrix
from .ztransform import is_multiplied_with
from copy import deepcopy
from warnings import warn


__all__ = ('DFT', 'DFTmatrix')


def is_in_interval(n0, n1, nk, comment=''):
    """Check if bin belongs to an interval [n0, n1].
     n0, n1, nk may contain parameter

    This returns:
     -1  if bin is not in interval and is on the left side
      0  if bin is in interval
      1  if bin is not in interval and is on the right side

    If n0, n1 or nk contain parameter, 0 is returned if no decision is
    possible and a comment is printed.

    """

    # Is in [n0, n1]
    if sym.sympify(nk - n0).is_nonnegative and sym.sympify(n1 - nk).is_nonnegative:
        return 0
    # Is on the left
    elif sym.sympify(nk - n0).is_negative:
        return -1
    # Is on the right
    elif sym.sympify(n1 - nk).is_negative:
        return 1
    # No decision, assume in interval
    else:
        warn('DFT of %s, assuming %s is in interval [%s, %s]' %
             (comment, nk, n0, n1))
        return 0


def find_pow(q_expr, base, NN):
    """This finds terms of the form base**NN in the expression q_expr.
    It is used for simplifying q**N terms to one if q=exp(-+j*2*pi/N*k).

    It returns a list of all bases found e.g., [q/3, q/5, 2q, aq,  .......].
    """

    result = []
    for x in q_expr.atoms(sym.Pow):
        this_ex = x.as_base_exp()[1]
        this_ba = x.as_base_exp()[0]
        if this_ex == NN and not sym.simplify(this_ba / base).has(base):
            result += [this_ba]

    return result


def special_diff(var, Q, o, s, var_a, var_N):
    """Make derivative of a**(k-s)*Q / (1-a**N)**o with respect to a, and
       divide result by a**(k-s-1)/(1-a**N)**(o+1) to obtain a
       polynom in n.  This is used with make_polynom_table.

    """
    return ((var - s) * Q + var_a*sym.diff(Q, var_a)) * (1 - var_a**var_N) + o * var_N * Q * var_a**var_N


def make_polynom_table1(var, MaxP, aa=None, NN=None, prn=False):
    """ This Table contains solutions of the DFT of
          exp(j*2*pi/N*n*mu) / (1-a*exp(j*2*pi/N*n))**p
          where mu < p  and mu, p are integer
          NO roots in the denominator

          use the table as follows
          the list Qpoly[p] contains polynoms with respect to k of the solutions to a specific p and  0<=mu<=p-1 .
          The full solution is obtained as
          N * a**(k-mu) / (1-a**N)**p / (p-1) *Qpoly[p][mu]           p >= 2
          N*a**k / (1-a**N)* TP[1]                                                       p =1

           Example for p=3
           1 / (1-a*exp(-j*2*pi/N*n))**3  = N*a**k * Qpoly[3][0] / (1-a**N)**3 / 2
           exp(-j*2*pi/N*n) / (1-a*exp(-j*2*pi/N*n))**3  = N*a**(k-1) * Qpoly[3][1] / (1-a**N)**3 / 2
           exp(-j*2*pi/N*n*2) / (1-a*exp(-j*2*pi/N*n))**3  = N*a**(k-2) * Qpoly[3][2] / (1-a**N)**3 / 2

          The results for all p's are calculated (--> Qpoly), only p=MaxP is returned.
          The parameter a and N can be set, if not a and N remain a symbol
    """
    a = sym.Symbol('a')
    N = sym.Symbol('N', integer=True)
    Qpoly = [[], [1]]
    # Make the Polynoms Q_p_mu
    # Run through all p, as the formula is recursiv (manuscript)
    for l in range(2, MaxP + 1, 1):
        # Special derivative when polynoms are involved (see also handle n**p in termXq)
        Q0 = special_diff(var, Qpoly[l - 1][0], l - 1, 0,
                          a, N) + Qpoly[l - 1][0] * (1 - a**N) * (l - 1)
        this_Qpoly = [sym.simplify(sym.expand(Q0))]
        # Run for each p through all 0<=mu<=p-1
        for m in range(l - 1):
            Qm = special_diff(var, Qpoly[l - 1][m], l - 1, m, a, N)
            this_Qpoly += [sym.simplify(sym.expand(Qm))]

        # Write result to a  list containing all polynomials
        Qpoly += [this_Qpoly]

    # Return polynoms to p=MaxP with specific a and N given
    ret_pol = []
    for mu in range(MaxP):
        pol = 0
        for r in range(MaxP):
            dj = sym.expand(Qpoly[MaxP][mu]).coeff(var, r)
            if not aa is None:
                dj = dj.subs(a, aa)
            if not NN is None:
                dj = dj.subs(N, NN)
            pol += dj * var**r
            if prn:
                print("mu = ", mu, " kn**", r, "  :  ", dj)
        if prn:
            print("Q = ", pol, "\n")
        ret_pol += [pol]
    # Return polynom for MaxP only, otherwise Qpoly must be returned
    return ret_pol


def make_polynom_table2(var, p, varN=None, prn=False):
    """This Table contains solutions of the DFT of
        exp(j*2*pi/N*n*mu) / (1-a*exp(j*2*pi/N*n))**p
        where mu < p  and mu, p are integer
        a is of the form exp(-j*2*pi*r/N), so here we have roots in the denominator at n=r

        The solution is of the form
            Q[mu] * exp(-j*2*pi*r/N*(k-mu)), were Q[mu] is a polynom in k of order p

            idea q = exp(j*2*pi/N*(n-r)):
           The transform of C_p(q) / (1-q)**p is (k**p - 1/N*sum_0_N-1(k**p)) * exp(-j*2*pi*mu*k/N)                 (1)

           Find from these the transforms of 1/(1-q)**p, q/(1-q)**p, ...., q**(p-1)/(1-q)**p by an appropriate
           linear combination of the transforms in (1).
           The C_p(q) is the polynom according to  sum_0_(N-1) n**p * q**n = N*Cp / (1-q)**p with q**N=1
    """
    NN = sym.Symbol('NN')
    q = sym.Symbol('q')
    a_p = 1 + 0 * q
    b_p = 1 + 0 * q
    AllC = [0]
    # Find first the C_p with the recursive formula for the repetitive derivative   q*diff((1-q**N)/(1-q), q)
    # by writing this as (a_p - q**N*b_p) /(1-q)**p.
    for i in range(p):
        a_p = sym.expand((q - q**2) * sym.diff(a_p, q) + (i + 1) * q * a_p)
        b_p = sym.expand(NN * (1 - q) * b_p + (q - q ** 2) *
                         sym.diff(b_p, q) + (i + 1) * q * b_p)
        pol, rem = sym.div(a_p - b_p, 1 - q)
        AllC += [(pol/NN).as_poly(q).as_expr()]
    if prn:
        print("p=1, ...., ", p, "  Ci(q) = ", AllC[1:])

    # Solve the equation x_p*C_p + x_(p-1)*C_(p-1)(1-q) +... = q**s = r_1 +r_2*q + .. r_p*q**(p-1)
    # to find the linear factors xi by equating coefficients in the variable q
    x = [sym.symbols('x%d' % i) for i in range(p)]
    r = [sym.symbols('r%d' % i) for i in range(p)]
    pol = 0 * q
    for i in range(p):
        # Left side of equation
        pol += x[i] * AllC[i + 1] * (1 - q) ** (p - i - 1)

    # Expand left side in orders of q
    coeffs = sym.Poly(sym.expand(pol), q).all_coeffs()[-1::-1]
    # print("Left side equation ", coeffs)
    eq = []
    xi = []
    # Subtract right hand side
    for i in range(p):
        eq += [coeffs[i] - r[i]]
        xi += [x[i]]

    # Solve for the xi by equating coefficients in q
    sol = sym.linsolve(eq, tuple(xi))
    All_rhs = sym.eye(p)
    # Set r_1, .., r_p to obtain the solution for all desired new transforms
    # e.g.  if q**2/(1-q)**p is to be found, then chose r_1=r_2=0    r_3=1   r_4=..=r_p=0
    ret_pol = []
    for i in range(p):
        rhs = list(All_rhs[:, i])
        repl = dict(zip(r, rhs))

        if prn:
            print(
                "q**", i, " :  linear factors of basic transform [x_i] = ", sol.subs(repl).args[0])
        # Find finally the polynoms in k for the new transforms  1/(1-q), , , , q**(p-1)/(1-q)**p
        trans_i = 0 * var
        for ii in range(p):
            # Sum k**(ii+1)
            faulh = sym.factor(
                (sym.bernoulli(ii + 2, NN) - sym.bernoulli(ii + 2, 0)) / (ii + 2))
            trans_i += sol.subs(repl).args[0][ii] * \
                (var ** (ii + 1) - faulh / NN)
        if not varN is None:
            trans_i = trans_i.subs(NN, varN)
        ret_pol += [trans_i.as_poly(var).as_expr()]
        if prn:
            print("q**", i, " :  D(.) = ", trans_i.as_poly(var).as_expr(), "\n")

    return ret_pol


# Does not use simplify as sometimes weird expression may be generated
# e.g. exp(j*2*pi/N) = (-1)**(2/N), or denoms are expanded and not
# factorized or it takes to long
def simp_rat(expr, q, method=1):
    """
    Make a few simplification steps for a rational function, useful for DFT's
    method 1: uses cancel
    method 2: uses factor_list of num and denom and cancel their common factors by hand
    """
    if method == 1:
        #ret_expr = sym.simplify(expr)
        #ret_expr = sym.cancel(expr, q)
        try:
            # that does not work with parameter
            ret_expr = sym.cancel(expr, q)
        except:
            ret_expr = sym.cancel(expr)
            return ret_expr
        Denom = sym.denom(ret_expr)
        # Denominator was probably expanded, so return original one
        if Denom.is_Add and len(Denom.args) == 3:
            return sym.numer(expr).as_poly(q).as_expr() / sym.denom(expr)
        else:
            ret_expr = sym.numer(ret_expr).as_poly(
                q).as_expr() / sym.factor(Denom, q)
            return ret_expr
    elif method == 2:
        Num = sym.numer(expr)
        Denom = sym.denom(expr)
        # Make factor list
        f_n = list(sym.factor_list(Num, q))
        f_d = list(sym.factor_list(Denom, q))
        # Check for common factors
        for ii, xx in enumerate(f_d[1]):
            for jj, yy in enumerate(f_n[1]):
                hits = sym.cancel(yy[0] / xx[0])
                if not hits.has(q):
                    # Common factor found
                    f_n[0] *= hits**yy[1]
                    # Exponents
                    p_n = yy[1]
                    p_d = xx[1]
                    # Cancel common factors and write as positive exponent
                    if p_n <= p_d:
                        f_n[1][jj] = (xx[0], 0)
                        f_d[1][ii] = (xx[0], p_d-p_n)
                    else:
                        f_n[1][jj] = (xx[0], p_n-p_d)
                        f_d[1][ii] = (xx[0], 0)
                    break

        # Rebuild factorization
        Denom = 1 + 0 * q
        Num = f_n[0] / f_d[0]
        for xx in f_d[1]:
            Denom *= xx[0]**xx[1]
        for xx in f_n[1]:
            Num *= xx[0]**xx[1]
        ret_expr = Num/Denom
        return ret_expr
    else:
        raise ValueError('Invalid method')


class QkTransform(object):
    """This handles X(q) transforms, where

    X(q) = sum_lower^upper x(n) * q**N

    As special cases are possible (e.g., for q=0) the class handles the
    general expression and special cases"""

    def __init__(self, val_q, k0=None, val0=None, valqk=None, val_k=0):
        self.Xq = val_q
        self.Xk = val_k
        if k0 is not None:
            self.has_special = True
            self.cases = [k0]
            self.case_expr = {k0: [val0, valqk]}
        else:
            self.has_special = False
            self.cases = []
            self.case_expr = {}

    def multiply(self, c0):
        """Multiply by constant."""

        self.Xq *= c0
        self.Xk *= c0
        if self.has_special:
            for key in self.case_expr:
                self.case_expr[key][0] *= c0
                self.case_expr[key][1] *= c0

    def rm_cases(self):
        """Delete special cases."""

        if self.has_special:
            self.has_special = False
            self.cases = []
            self.case_expr = {}

    def shift_k(self, dk):
        """Shift the special cases k values by dk."""

        if self.has_special:
            self.cases = [ii + dk for ii in self.cases]
            new_dict = {}
            for key in self.case_expr:
                new_dict[key + dk] = self.case_expr[key]
            self.case_expr = new_dict

    def add(self, sec):
        """Merge (add) two QkTransforms."""

        self.Xq += sec.Xq
        self.Xk += sec.Xk
        # Both have special cases
        if sec.has_special and self.has_special:
            self.cases += [ii for ii in sec.cases if ii not in self.cases]
            for key in sec.case_expr:
                if key in self.case_expr:
                    self.case_expr[key][0] += sec.case_expr[key][0]
                    self.case_expr[key][1] += sec.case_expr[key][1]
                else:
                    self.case_expr[key] = sec.case_expr[key]
        # Only the second has special cases
        elif sec.has_special and not self.has_special:
            self.has_special = True
            self.cases = [ii for ii in sec.cases]
            for key in sec.case_expr:
                self.case_expr[key] = sec.case_expr[key]

    def subs(self, q, q_new):
        """Substitute q by a new expression."""

        self.Xq = (self.Xq).subs(q, q_new)
        if self.has_special:
            for key in self.case_expr:
                self.case_expr[key][1] = self.case_expr[key][1].subs(q, q_new)

    def simp_qN(self, ba, NN, all=False):
        """Simplifies ba**NN to 1

        If all is False, only the special cases will regarded.
        Replace e.g., (q/3)**NN by (1/3)**N
        """

        # Simplify for "no special cases" or if demanded
        if all or not self.has_special:
            # General expresion
            terms = find_pow(self.Xq, ba, NN)
            for full_base in terms:
                self.Xq = self.Xq.replace(full_base**NN, (full_base / ba)**NN)

        # Special cases
        if self.has_special:
            for key in self.case_expr:
                terms = find_pow(self.case_expr[key][1], ba, NN)
                for full_base in terms:
                    self.case_expr[key][1] = self.case_expr[key][1].replace(
                        full_base**NN, (full_base / ba)**NN)

    def make_transform(self, NN, backward, q, sub, k, piecewise):
        """Assemble the final tranform from all the cases, the general case, or Xk.

        backward is True for IDFT
        q is replaced by sub
        """
        simp_method = 1  # use cancel
        # simp_method = 2  # use home made factorization with cancellation
        if self.has_special:
            for ca in self.cases:
                k0 = ca
                # Shift k0 to [0, N-1]
                if backward and sym.sympify(ca).is_positive:
                    k0 -= NN
                elif not backward and sym.sympify(ca).is_negative:
                    k0 += NN

                arg_delta = k - k0
                if backward:
                    # as k=-n and k0 is negative
                    arg_delta *= -1

                # Special case * delta
                self.Xk += self.case_expr[ca][0] * UnitImpulse(arg_delta)

                # Rest * (1-delta) to avoid unwanted superposition try
                # to make expressions as compact as possible by using
                # simp_rat
                this_expr = self.case_expr[ca][1]
                if this_expr.is_Add:
                    sum_expr = 0 * q
                    for expr_i in this_expr.args:
                        sum_expr += simp_rat(expr_i, q, method=simp_method)
                    this_expr = sum_expr
                else:
                    this_expr = simp_rat(this_expr, q, method=simp_method)
                this_expr = this_expr.subs(q, sub)

                # Check if val is at k0 is unlimted and prepare for piecewise
                if not backward:
                    case_val = this_expr.subs(k, k0)
                    neq = sym.Ne(k, k0)
                else:
                    case_val = this_expr.subs(-k, -k0)
                    neq = sym.Ne(-k, -k0)

                if (case_val == sym.nan or case_val.has(sym.zoo) or
                        case_val == sym.oo or case_val == -sym.oo):

                    if piecewise:
                        self.Xk += sym.Piecewise((this_expr, neq))
                    else:
                        self.Xk += this_expr * (1 - UnitImpulse(arg_delta))
                else:
                    # Finite value, so write slightly differently
                    self.Xk += this_expr - case_val*UnitImpulse(arg_delta)
        else:
            this_expr = self.Xq
            if this_expr.is_Add:
                sum_expr = 0 * q
                for expr_i in this_expr.args:
                    sum_expr += simp_rat(expr_i, q, method=simp_method)
                this_expr = sum_expr
            else:
                this_expr = simp_rat(this_expr, q, method=simp_method)
            self.Xk = this_expr.subs(q, sub) + self.Xk


class DFTTransformer(BilateralForwardTransformer):

    name = 'DFT'
    is_inverse = False

    def key(self, expr, n, k, **kwargs):
        return expr, n, k, kwargs.get('N', None), kwargs.get('piecewise', False)

    def noevaluate(self, expr, n, k):

        foo = expr * sym.exp(-2 * j * pi * n * k / self.N)
        result = sym.Sum(foo, (n, 0, self.N - 1))
        return result

    def check(self, expr, n, k, N, **kwargs):

        try:
            N = N.expr
        except:
            pass

        # Convert constants to SymPy expressions.
        N = sym.sympify(N)

        if not N.is_integer and not N.is_positive:
            raise ValueError(
                "%s not positive integer, redefine as %s = symbol('%s', integer=True, positive=True)" % (N, N, N))

        symbols = expr.free_symbols
        symbol_names = [str(e) for e in symbols]
        if str(N) in symbol_names and N not in symbols:
            warn('There is a symbol in the expression with the same name as the DFT size %s but is not the same symbol' % N)

        self.N = N

        if expr.has(k):
            self.error('Expression depends on k')

        if expr.is_Piecewise and expr.args[0].args[1].has(n >= 0):
            self.error('Expression is unknown for n < 0 (use causal=True)')

    def sympy(self, expr, n, k):

        foo = expr * sym.exp(-2 * j * pi * n * k / self.N)
        result = sym.summation(foo, (n, 0, self.N - 1))

        return result

    def func(self, expr, n, k):

        if not isinstance(expr, AppliedUndef):
            self.error('Expecting function')

        scale, shift = scale_shift(expr.args[0], n)

        # Convert v(n) to V(k), etc.
        name = expr.func.__name__
        if self.is_inverse:
            func = sym.Function(name[0].lower() + name[1:])
        else:
            func = sym.Function(name[0].upper() + name[1:])

        result = func(k / scale) / abs(scale)

        if shift != 0:
            if self.is_inverse:
                shift = -shift
            result = result * \
                sym.exp(2 * sym.I * sym.pi * k * shift / (scale * self.N))

        if self.is_inverse:
            result *= self.N

        return result

    def function(self, expr, n, k):

        # Handle expressions with a function of FOO, e.g.,
        # v(n), v(n) * y(n), 3 * v(n) / n, v(4 * a * n), etc.,

        if not expr.has(AppliedUndef):
            self.error()

        const, expr = factor_const(expr, n)

        if isinstance(expr, AppliedUndef):
            return self.func(expr, n, k) * const

        tsym = sympify(str(n))
        expr = expr.subs(tsym, n)

        rest = sym.S.One
        undefs = []
        for factor in expr.as_ordered_factors():
            if isinstance(factor, AppliedUndef):
                if factor.args[0] != n:
                    self.error('Weird function %s not of %s' % (factor, n))
                undefs.append(factor)
            else:
                rest *= factor

        if rest.has(AppliedUndef):
            # Have something like 1/v(n)
            self.error()

        exprs = undefs
        if rest.has(n):
            exprs = exprs + [rest]
            rest = sym.S.One

        result = self.term(exprs[0], n, k) * rest

        if len(exprs) == 1:
            return result * const

        dummy = 'm' if self.is_inverse else 'l'

        for m in range(len(exprs) - 1):
            if m == 0:
                nu = miscsymbol(dummy, integer=True)
            else:
                nu = miscsymbol(dummy + '_%d' % m, integer=True)
            expr2 = self.term(exprs[m + 1], n, k)
            # Should be a circular convolution.
            result = sym.Sum(result.subs(k, k - nu) * expr2.subs(k, nu),
                             (nu, 0, self.N - 1)) / self.N

        return result * const

    # Make transform  Xq = sum_lower^upper  x[n] * q**n
    def termXq(self, expr, n, k, q, lower, upper):

        const, expr = factor_const(expr, n)
        args = expr.args
        xn_fac = []

        # Check for constant.
        if not expr.has(n):
            result_q = const * expr * q**lower * \
                (1 - q**(upper - lower + 1)) / (1 - q)
            # Special case k = 0
            result_1 = const * expr * (upper - lower + 1)
            result = QkTransform(result_q, 0, result_1, result_q)
            return result

        # Handle delta(n-n0)
        elif (expr.is_Function and expr.func == UnitImpulse and
              ((expr.args[0]).as_poly(n)).is_linear):
            aa = args[0].coeff(n, 1)
            bb = args[0].coeff(n, 0)
            nn0 = -bb / aa
            if nn0.is_integer:
                if is_in_interval(lower, upper, nn0, comment='delta') == 0:
                    # Shift frequency to -pi/2 ...  pi/2
                    if nn0.has(self.N):
                        nn0 = nn0.subs(self.N, 0)
                    elif sym.sympify(nn0 - self.N / 2).is_positive:
                        nn0 -= self.N

                    result_q = const * q**nn0
                    result = QkTransform(result_q)
                else:
                    result = QkTransform(0 * q)
                return result
            else:
                print("delta: bin ", nn0, " is not an integer")

        # Handle n**p
        if (expr == n or
            (expr.is_Pow and args[1].is_integer and args[1].is_positive
             and args[0] == n)):
            p = 1
            try:
                p = args[1]
            except:
                pass
            # Derivatives
            result_q = (q**lower - q**(upper + 1)) / (1 - q)
            # Make the first few by hand to make expressions shorter
            # and faster, use in all cases for the derivatives (and
            # multiplication with q) the structure
            # (q**lower * polynomial (q)- q**(upper+1)*polynomial(q)) (1-q)**(p+1) . The order of the polynomials is equal to p
            if p == 1:
                A_l = lower + q * (1 - lower)
                B_u = -q * upper + upper + 1
            elif p == 2:
                A_l = lower**2 + q**2 * \
                    (lower**2 - 2 * lower + 1) + q * \
                    (-2 * lower**2 + 2 * lower + 1)
                B_u = q**2*upper**2 + q * \
                    (-2 * upper**2 - 2 * upper + 1) + upper**2 + 2 * upper + 1
            elif p == 3:
                A_l = lower**3 + q**3 * (-lower**3 + 3 * lower**2 - 3 * lower + 1) + q**2*(
                    3 * lower**3 - 6*lower**2 + 4) + q*(-3 * lower**3 + 3 * lower**2 + 3 * lower + 1)
                B_u = -q**3 * upper**3 + q**2*(3 * upper**3 + 3 * upper**2 - 3 * upper + 1) + q*(
                    -3 * upper**3 - 6*upper**2 + 4) + upper**3 + 3 * upper**2 + 3 * upper + 1
            # Generate the polynoms A_l, B_u individually
            # TODO additional option to make a smart simplification by inspection of upper and lower
            else:
                # First polynoms
                A_l = 1 + 0 * q
                B_u = 1 + 0 * q
                # Find the next recursively
                for i in range(p):
                    A_l = sym.expand(
                        lower * A_l + (q - q**2) * sym.diff(A_l, q) - q * lower * A_l + (i + 1) * q * A_l)
                    B_u = sym.expand((upper + 1) * B_u + (q - q**2) * sym.diff(
                        B_u, q) - (upper + 1) * B_u * q + (i + 1) * q * B_u)
                A_l = A_l.as_poly(q).as_expr()
                B_u = B_u.as_poly(q).as_expr()
            result_q = const * (q**lower * A_l - q **
                                (upper + 1) * B_u) / (1 - q)**(p + 1)

            # Special case k=0, use Faulhaber's formula
            result_1 = const * \
                sym.factor((sym.bernoulli(p + 1, upper + 1) -
                           sym.bernoulli(p + 1, lower)) / (p + 1))
            result = QkTransform(result_q, 0, result_1, result_q)
            return result

        # Handle * rect((n-a)/b)
        elif is_multiplied_with(expr, n, 'rect', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = 1 / sym.expand(ref[0]).coeff(n, 1)
            aa = -bb * sym.expand(ref[0]).coeff(n, 0)
            # Left and right index
            nn0 = aa - bb // 2
            nn1 = nn0 + bb - 1
            if (not aa.is_integer) or (not bb.is_integer):
                print("rect((n-a)/b): parameter not an integer")
            else:
                r1 = is_in_interval(lower, upper, nn0,
                                    comment='rect left step')
                r2 = is_in_interval(lower, upper, nn1,
                                    comment='rect right step')
                if r1 == 0 and r2 == 0:
                    result = self.termXq(expr, n, k, q, nn0, nn1)
                elif r1 == -1 and r2 == 0:
                    result = self.termXq(expr, n, k, q, lower, nn1)
                elif r1 == 0 and r2 == 1:
                    result = self.termXq(expr, n, k, q, nn0, upper)
                elif r1 == -1 and r2 == 1:
                    result = self.termXq(expr, n, k, q, lower, upper)
                else:
                    result = QkTransform(0 * q)

                result.multiply(const)
                return result

        # Handle * u(n-n0)
        elif is_multiplied_with(expr, n, 'UnitStep', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            aa = ref[0].coeff(n, 1)
            bb = ref[0].coeff(n, 0)
            if abs(aa) != 1:
                print("Use u(n-n0)")
            if not bb.is_integer:
                print("Step(n-n0): n0 not an integer")
            else:
                # Positive step
                if aa == 1:
                    nn0 = -bb
                    rg = is_in_interval(lower, upper, nn0, comment='step')
                    if rg == -1:
                        result = self.termXq(expr, n, k, q, lower, upper)
                    elif rg == 0:
                        result = self.termXq(expr, n, k, q, nn0, upper)
                    else:
                        result = QkTransform(0 * q)

                # Negative step
                elif aa == -1:
                    nn0 = bb
                    rg = is_in_interval(lower, upper, nn0, comment='Step:')
                    if rg == 1:
                        result = self.termXq(expr, n, k, q, lower, upper)
                    elif rg == 0:
                        result = self.termXq(expr, n, k, q, lower, nn0)
                    else:
                        result = QkTransform(0 * q)

                result.multiply(const)
                return result

        # Handle * exp(j*a*n+b)
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac) and abs(xn_fac[-1] / sym.exp(args[0].coeff(n, 0))) == 1:
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            aa = sym.expand(ref[0]).coeff(n, 1) / sym.I
            bb = sym.expand(ref[0]).coeff(n, 0)
            # Find transform
            result = self.termXq(expr, n, k, q, lower, upper)
            # Check frequency
            if aa.is_constant() and abs(aa) > pi:
                warn("Frequency may be out of range")

            result.subs(q, q * sym.exp(sym.I * aa))
            # Check special case and shift accordingly
            k0 = aa * self.N / 2 / pi
            if k0.is_integer and result.has_special:
                result.shift_k(k0)
            else:
                result.rm_cases()
            result.multiply(const * sym.exp(bb))
            return result

        # Handle * sin(b*n+c)
        elif is_multiplied_with(expr, n, 'sin(n)', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)
            # Check frequency
            if bb.is_constant() and abs(bb) > pi:
                warn("Frequency may be out of range")

            result = self.termXq(expr, n, k, q, lower, upper)
            # Make copies for the transformation:
            # X(q) = (exp(j*bb)*Xq (q*exp(j*bb)) - exp(-j*bb)*Xq(q*exp(-j*bb))) / 2 / j
            rq1 = deepcopy(result)
            rq2 = deepcopy(result)
            # Make general shift
            rq1.subs(q, q * sym.exp(sym.I * bb))
            rq2.subs(q, q * sym.exp(-sym.I * bb))
            # Check special case shift s
            k0 = bb * self.N / 2 / pi
            if k0.is_integer and result.has_special:
                rq1.shift_k(k0)
                rq2.shift_k(-k0)
            else:
                rq1.rm_cases()
                rq2.rm_cases()

            rq1.multiply(const * sym.exp(sym.I * cc) / 2 / sym.I)
            rq2.multiply(-const * sym.exp(-sym.I * cc) / 2 / sym.I)

            # Add both parts
            rq1.add(rq2)
            return rq1

        # Handle * cos(b*n+c)
        elif is_multiplied_with(expr, n, 'cos(n)', xn_fac):
            expr /= xn_fac[-1]
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)
            # Check frequency
            if bb.is_constant() and abs(bb) > pi:
                warn("Frequency may be out of range")

            result = self.termXq(expr, n, k, q, lower, upper)
            # Make copies for the transformation:
            # X(q) = (exp(j*bb)*Xq (q*exp(j*bb)) + exp(-j*bb)*Xq(q*exp(-j*bb))) / 2
            rq1 = deepcopy(result)
            rq2 = deepcopy(result)
            # Make general shifts
            rq1.subs(q, q * sym.exp(sym.I * bb))
            rq2.subs(q, q * sym.exp(-sym.I * bb))

            # Check special case shifts
            k0 = bb * self.N / 2 / pi
            if k0.is_integer and result.has_special:
                rq1.shift_k(k0)
                rq2.shift_k(-k0)
            else:
                rq1.rm_cases()
                rq2.rm_cases()
            rq1.multiply(const * sym.exp(sym.I * cc) / 2)
            rq2.multiply(const * sym.exp(-sym.I * cc) / 2)

            # Add both parts
            rq1.add(rq2)
            return rq1

        # Handle * exp(bb*n+cc)
        elif is_multiplied_with(expr, n, 'exp(n)', xn_fac):
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            bb = ref[0].coeff(n, 1)
            cc = ref[0].coeff(n, 0)
            result = self.termXq(expr, n, k, q, lower, upper)
            result.subs(q, q * sym.exp(bb))
            # No special cases remain
            result.rm_cases()
            result.multiply(const * sym.exp(cc))
            return result

        # Handle * a**n
        elif is_multiplied_with(expr, n, 'a**n', xn_fac):
            expr /= xn_fac[-1]
            expr = sym.simplify(expr)
            ref = xn_fac[-1].args
            lam = ref[0]
            bb = ref[1].coeff(n, 1)
            cc = ref[1].coeff(n, 0)
            result = self.termXq(expr, n, k, q, lower, upper)
            result.subs(q, q * lam**bb)
            # No special cases remain
            result.rm_cases()
            result.multiply(const * lam**cc)
            return result

        # Handle * n
        elif is_multiplied_with(expr, n, 'n', xn_fac):
            expr = expr / xn_fac[-1]
            result = self.termXq(expr, n, k, q, lower, upper)
            result.Xq = q * sym.diff(result.Xq, q)
            if result.has_special:
                raise ValueError(
                    "No sym.diff possible for discrete cases, refine handles")
            else:
                return result

        # No case matches so return None as this can be used to check if
        # a case was found and if a simplify is useful for the general
        # summation term simplify runs too long
        return None

    def termXk(self, expr, n, k):
        """Make transform  X = sum_0^N-1  x[n] * exp(-2*pi*j*k*n/N)
           BUT result can not be written as function of exp(-2*pi*j*k/N)
           This part will be mainly important and used with the inverse DFT

           In general the terms implemented have the form:

              exp(j*2*pi/N*n*mu)/(1-a*exp(j*2*pi/N*n))**p      with mu integer

              or exp(-j*2*pi/N*k*mu)/(1-a*exp(-j*2*pi/N*k))**p
        """

        const, expr = factor_const(expr, n)

        qq = sym.Symbol('qq')
        p0 = sym.Wild('p0', properties=[
                      lambda r: r.is_integer and r.is_nonnegative])

        if not self.is_inverse:
            sg = 1
        else:
            sg = -1

        # Replace exp(j*2*pi/N*n) by qq (forward) on calling termXk
        # for the inverse, n is k and k is n
        expr_qq = expr.replace(
            sym.exp(sg*sym.I * 2 * pi * n / self.N * p0), qq**p0)
        # kn is (always) the output variable
        kn = k

        # TODO make polynom in qq, 1/qq terms may occur?

        # First case, Denom is never zero (i.e. for no integer n)
        if not expr_qq.has(n) and expr_qq.is_rational_function(qq):
            # Find numerator and denominator
            Num = sym.numer(expr_qq)
            Denom = sym.denom(expr_qq)

            # Polynomial division  d0+d1*q+..+dN*q**N + rem/Denom
            pol, rem = sym.div(Num, Denom)
            poly_coeff = sym.Poly(pol, qq).all_coeffs()[-1::-1]

            # First part c0+c1*q+....
            result_0 = 0
            # Polynom
            for i, xx in enumerate(poly_coeff):
                result_0 += self.N * xx * UnitImpulse(kn-i) / 2

            # Second part
            #   Denominator is linear D=1-a*qq
            if Denom.as_poly(qq).is_linear:
                a0 = sym.expand(Denom).coeff(qq, 1)
                c0 = sym.expand(Denom).coeff(qq, 0)
                # Remainder
                aa = -a0 / c0
                result = rem * self.N * aa**k / (1 - aa**self.N) / c0 / 2
                return 2 * const * (result_0 + result)

            # Denominator is D=(1-a*qq)**p, p>=2
            elif Denom.is_Pow and Denom.args[0].as_poly(qq).is_linear and Denom.args[1].is_integer and Denom.args[1] >= 2:
                a0 = sym.expand(Denom.args[0]).coeff(qq, 1)
                c0 = sym.expand(Denom.args[0]).coeff(qq, 0)
                aa = -a0 / c0

                aN = aa**self.N
                NN = self.N
                pp = Denom.args[1]
                # The polynomials in k of order pp-1 of the solution (see manuscript)
                if pp == 2:
                    Q_p_mu = [(NN-1)*aN + 1 + kn*(1 - aN), NN*aN + kn*(1 - aN)]
                elif pp == 3:
                    a2N = aN**2
                    N2 = NN**2
                    Q_p_mu = [N2*a2N + N2*aN - 3*NN*a2N + 3*NN*aN + 2*a2N - 4*aN + 2 + kn * (-2*NN*a2N + 2*NN*aN + 3*a2N - 6*aN + 3) + kn**2*(aN - 1)**2,
                              N2*a2N + N2*aN - NN*a2N + NN*aN + kn *
                              (-2*NN*a2N + 2*NN*aN + a2N -
                               2*aN + 1) + kn**2*(aN - 1)**2,
                              N2*a2N + N2*aN + NN*a2N - NN*aN + kn*(-2*NN*a2N + 2*NN*aN - a2N + 2*aN - 1) + kn**2*(aN-1)**2]
                elif pp == 4:
                    a2N = aN**2
                    a3N = aN**3
                    N2 = NN**2
                    N3 = NN**3
                    Q_p_mu = [N3*a3N + 4*N3*a2N + N3*aN - 6*N2*a3N + 6*N2*aN + 11*NN*a3N - 22*NN*a2N + 11*NN*aN - 6*a3N + 18*a2N - 18*aN + 6 +
                              kn**3*(-a3N + 3*a2N - 3*aN + 1) +
                              kn**2*(3*NN*a3N - 6*NN*a2N + 3*NN*aN - 6*a3N + 18*a2N - 18*aN + 6) +
                              kn*(-3*N2*a3N + 3*N2*aN + 12*NN*a3N - 24*NN *
                                  a2N + 12*NN*aN - 11*a3N + 33*a2N - 33*aN + 11),
                              N3*a3N + 4*N3*a2N + N3*aN - 3*N2*a3N + 3*N2*aN + 2*NN*a3N - 4*NN*a2N + 2*NN*aN +
                              k**3*(-a3N + 3*a2N - 3*aN + 1) +
                              k**2*(3*NN*a3N - 6*NN*a2N + 3*NN*aN - 3*a3N + 9*a2N - 9*aN + 3) +
                              k*(-3*N2*a3N + 3*N2*aN + 6*NN*a3N - 12*NN *
                                 a2N + 6*NN*aN - 2*a3N + 6*a2N - 6*aN + 2),
                              N3*a3N + 4*N3*a2N + N3*aN - NN*a3N + 2*NN*a2N - NN*aN +
                              k**3*(-a3N + 3*a2N - 3*aN + 1) +
                              k**2*(3*NN*a3N - 6*NN*a2N + 3*NN*aN) +
                              k*(-3*N2*a3N + 3*N2*aN + a3N - 3*a2N + 3*aN - 1),
                              N3*a3N + 4*N3*a2N + N3*aN + 3*N2*a3N - 3*N2*aN + 2*NN*a3N - 4*NN*a2N + 2*NN*aN +
                              k**3*(-a3N + 3*a2N - 3*aN + 1) +
                              k**2*(3*NN*a3N - 6*NN*a2N + 3*NN*aN + 3*a3N - 9*a2N + 9*aN - 3) +
                              k*(-3*N2*a3N + 3*N2*aN - 6*NN*a3N + 12*NN*a2N - 6*NN*aN - 2*a3N + 6*a2N - 6*aN + 2)]
                elif pp == 5:
                    a2N = aN**2
                    a3N = aN**3
                    a4N = aN**4
                    N2 = NN**2
                    N3 = NN**3
                    N4 = NN**4
                    Q_p_mu = [N4*a4N + 11*N4*a3N + 11*N4*a2N + N4*aN - 10*N3*a4N - 30*N3*a3N + 30*N3*a2N + 10*N3*aN + 35*N2*a4N - 35*N2*a3N - 35*N2*a2N + 35*N2*aN - 50*NN*a4N + 150*NN*a3N - 150*NN*a2N + 50*NN*aN + 24*a4N - 96*a3N + 144*a2N - 96*aN + 24 +
                              kn**4*(a4N - 4*a3N + 6*a2N - 4*aN + 1) +
                              kn**3*(-4*NN*a4N + 12*NN*a3N - 12*NN*a2N + 4*NN*aN + 10*a4N - 40*a3N + 60*a2N - 40*aN + 10) +
                              kn**2*(6*N2*a4N - 6*N2*a3N - 6*N2*a2N + 6*N2*aN - 30*NN*a4N + 90*NN*a3N - 90*NN*a2N + 30*NN*aN + 35*a4N - 140*a3N + 210*a2N - 140*aN + 35) +
                              kn*(-4*N3*a4N - 12*N3*a3N + 12*N3*a2N + 4*N3*aN + 30*N2*a4N - 30*N2*a3N - 30*N2*a2N + 30*N2 *
                                  aN - 70*NN*a4N + 210*NN*a3N - 210*NN*a2N + 70*NN*aN + 50*a4N - 200*a3N + 300*a2N - 200*aN + 50),
                              N4*a4N + 11*N4*a3N + 11*N4*a2N + N4*aN - 6*N3*a4N - 18*N3*a3N + 18*N3*a2N + 6*N3*aN + 11*N2*a4N - 11*N2*a3N - 11*N2*a2N + 11*N2*aN - 6*NN*a4N + 18*NN*a3N - 18*NN*a2N + 6*NN*aN + k**4*(a4N - 4*a3N + 6*a2N - 4*aN + 1) +
                              kn**3*(-4*NN*a4N + 12*NN*a3N - 12*NN*a2N + 4*NN*aN + 6*a4N - 24*a3N + 36*a2N - 24*aN + 6) +
                              kn**2*(6*N2*a4N - 6*N2*a3N - 6*N2*a2N + 6*N2*aN - 18*NN*a4N + 54*NN*a3N - 54*NN*a2N + 18*NN*aN + 11*a4N - 44*a3N + 66*a2N - 44*aN + 11) +
                              kn*(-4*N3*a4N - 12*N3*a3N + 12*N3*a2N + 4*N3*aN + 18*N2*a4N - 18*N2*a3N - 18*N2*a2N + 18 *
                                  N2*aN - 22*NN*a4N + 66*NN*a3N - 66*NN*a2N + 22*NN*aN + 6*a4N - 24*a3N + 36*a2N - 24*aN + 6),
                              N4*a4N + 11*N4*a3N + 11*N4*a2N + N4*aN - 2*N3*a4N - 6*N3*a3N + 6*N3*a2N + 2*N3*aN - N2*a4N + N2*a3N + N2*a2N - N2*aN + 2*NN*a4N - 6*NN*a3N + 6*NN*a2N - 2*NN*aN +
                              kn**4*(a4N - 4*a3N + 6*a2N - 4*aN + 1) +
                              kn**3*(-4*NN*a4N + 12*NN*a3N - 12*NN*a2N + 4*NN*aN + 2*a4N - 8*a3N + 12*a2N - 8*aN + 2) +
                              kn**2*(6*N2*a4N - 6*N2*a3N - 6*N2*a2N + 6*N2*aN - 6*NN*a4N + 18*NN*a3N - 18*NN*a2N + 6*NN*aN - a4N + 4*a3N - 6*a2N + 4*aN - 1) +
                              kn*(-4*N3*a4N - 12*N3*a3N + 12*N3*a2N + 4*N3*aN + 6*N2*a4N - 6*N2*a3N - 6*N2*a2N + 6 *
                                  N2*aN + 2*NN*a4N - 6*NN*a3N + 6*NN*a2N - 2*NN*aN - 2*a4N + 8*a3N - 12*a2N + 8*aN - 2),
                              N4*a4N + 11*N4*a3N + 11*N4*a2N + N4*aN + 2*N3*a4N + 6*N3*a3N - 6*N3*a2N - 2*N3*aN - N2*a4N + N2*a3N + N2*a2N - N2*aN - 2*NN*a4N + 6*NN*a3N - 6*NN*a2N + 2*NN*aN +
                              kn**4*(a4N - 4*a3N + 6*a2N - 4*aN + 1) + k**3*(-4*NN*a4N + 12*NN*a3N - 12*NN*a2N + 4*NN*aN - 2*a4N + 8*a3N - 12*a2N + 8*aN - 2) +
                              kn**2*(6*N2*a4N - 6*N2*a3N - 6*N2*a2N + 6*N2*aN + 6*NN*a4N - 18*NN*a3N + 18*NN*a2N - 6*NN*aN - a4N + 4*a3N - 6*a2N + 4*aN - 1) +
                              kn*(-4*N3*a4N - 12*N3*a3N + 12*N3*a2N + 4*N3*aN - 6*N2*a4N + 6*N2*a3N + 6*N2*a2N - 6 *
                                  N2*aN + 2*NN*a4N - 6*NN*a3N + 6*NN*a2N - 2*NN*aN + 2*a4N - 8*a3N + 12*a2N - 8*aN + 2),
                              N4*a4N + 11*N4*a3N + 11*N4*a2N + N4*aN + 6*N3*a4N + 18*N3*a3N - 18*N3*a2N - 6*N3*aN + 11*N2*a4N - 11*N2*a3N - 11*N2*a2N + 11*N2*aN + 6*NN*a4N - 18*NN*a3N + 18*NN*a2N - 6*NN*aN +
                              kn**4*(a4N - 4*a3N + 6*a2N - 4*aN + 1) + k**3*(-4*NN*a4N + 12*NN*a3N - 12*NN*a2N + 4*NN*aN - 6*a4N + 24*a3N - 36*a2N + 24*aN - 6) +
                              kn**2*(6*N2*a4N - 6*N2*a3N - 6*N2*a2N + 6*N2*aN + 18*NN*a4N - 54*NN*a3N + 54*NN*a2N - 18*NN*aN + 11*a4N - 44*a3N + 66*a2N - 44*aN + 11) +
                              kn*(-4*N3*a4N - 12*N3*a3N + 12*N3*a2N + 4*N3*aN - 18*N2*a4N + 18*N2*a3N + 18*N2*a2N - 18*N2*aN - 22*NN*a4N + 66*NN*a3N - 66*NN*a2N + 22*NN*aN - 6*a4N + 24*a3N - 36*a2N + 24*aN - 6)]

                else:
                    # TODO  make faster
                    # we increase performance for pp >= 6 by saving the table in an extra file and read the polynoms and subs
                    # or make make_polynom_table1 faster by substituting a**(r*N) by arN and so on
                    Q_p_mu = make_polynom_table1(kn, pp, aa, NN)

                # case   d0/(1-a*qq)**p   + d1*qq/(1-a*qq)**p + .... + dm *qq**(p-1)/(1-a*qq)**p
                num_coeff = sym.Poly(rem, qq).all_coeffs()[-1::-1]
                if len(num_coeff) > len(Q_p_mu):
                    raise ValueError(
                        " Error : polynomial division may be wrong! ")
                # sum up
                result = 0
                for i in range(len(num_coeff)):
                    result += num_coeff[i] * Q_p_mu[i] / aa**i
                result *= NN / \
                    sym.factorial(pp - 1) * aa**kn / \
                    (1 - aN) ** pp / c0**pp / 2
                return 2 * const * (result_0 + result)

            else:
                return None

        # Second case, denominators has roots for integer n,
        # so term must contain a factor (1-delta(k-k0)) to indicate that k!=k0
        # better the expression is defined with sym.Piecewise sym.Piecewise((expr, sym.Ne(k, k0)))
        # search for a factor (1-delta(n-n0)) or for piecewise.
        this_case = False
        # Check for piecewise
        if expr.is_Piecewise:
            match = list(expr.args[0][1].find(sym.Ne(n, p0)))
            n0 = match[0].args[1]
            div_expr = expr_qq.args[0][0]
            this_case = True
        # Check for (1-delta) factor
        else:
            match = list(expr.find(1 - UnitImpulse(n - p0)))
            if len(match) == 1:
                n0 = -((-match[0] + 1).args[0] - n)
                div_expr = expr_qq / match[0]
                this_case = True
        if this_case:
            if not div_expr.has(n) and div_expr.is_rational_function(qq):
                # Find numerator and denominator
                Num = sym.numer(div_expr)
                Denom = sym.denom(div_expr)

                # Polynomial division  d0+d1*q+..+dN*q**N + rem/N
                pol, rem = sym.div(Num, Denom)
                poly_coeff = sym.Poly(pol, qq).all_coeffs()[-1::-1]

                # First part c0+c1*q+....
                result_0 = 0
                # Polynom
                for i, xx in enumerate(poly_coeff):
                    # as n!=n0 this n has to be excluded
                    if i != n0:
                        result_0 += xx * UnitImpulse(kn - i) * self.N

                # Denominator is linear D=1-a*qq
                if Denom.as_poly(qq).is_linear:
                    a0 = sym.expand(Denom).coeff(qq, 1)
                    c0 = sym.expand(Denom).coeff(qq, 0)
                    # remainder
                    aa = -a0 / c0
                    if not aa == sym.exp(-sg * j * 2 * pi / self.N * n0):
                        print("Rewrite expression as ", expr /
                              match[0] - expr.subs(n, n0))
                        return None
                    # Result
                    result = ((self.N - 1) / 2 - kn) * \
                        sym.exp(-sg * j * 2 * pi * k * n0 / self.N) / c0
                    return const * (result_0 + result)

                # Denominator is D=(1-a*qq)**p, p>=2
                elif Denom.is_Pow and Denom.args[0].as_poly(qq).is_linear and Denom.args[1].is_integer and Denom.args[1] >= 2:
                    a0 = sym.expand(Denom.args[0]).coeff(qq, 1)
                    c0 = sym.expand(Denom.args[0]).coeff(qq, 0)
                    aa = -a0 / c0
                    if not aa == sym.exp(-sg * j * 2 * pi / self.N * n0):
                        print("Rewrite expression as ", expr /
                              match[0] - expr.subs(n, n0))
                        return None

                    NN = self.N
                    pp = Denom.args[1]
                    # The polynomials in k of order pp of the solution (see manuscript)
                    if pp == 2:
                        T_p_mu = [(-NN**2 + 6*NN - 6*kn**2 + 6*kn*(NN - 2) - 5) / 12,
                                  (-NN**2 + 6*NN*kn - 6*kn**2 + 1) / 12]
                    elif pp == 3:
                        T_p_mu = [(-3*NN**2 + 12*NN - 4*kn**3 + 6*kn**2*(NN - 3) + kn*(-2*NN**2 + 18*NN - 24) - 9) / 24,
                                  (-NN**2 - 4*kn**3 + 6*kn**2*(NN - 1) +
                                   kn*(-2*NN**2 + 6*NN) + 1) / 24,
                                  (NN**2 - 4*kn**3 + 6*kn**2*(NN + 1) + kn*(-2*NN**2 - 6*NN) - 1) / 24]
                    elif pp == 4:
                        T_p_mu = [(NN**4 - 110*NN**2 + 360*NN - 30*kn**4 + kn**3*(60*NN - 240) + kn**2*(-30*NN**2 + 360*NN - 660) + kn*(-120*NN**2 + 660*NN - 720) - 251) / 720,
                                  (NN**4 - 20*NN**2 - 30*kn**4 + kn**3*(60*NN - 120) + kn**2 *
                                   (-30*NN**2 + 180*NN - 120) + kn*(-60*NN**2 + 120*NN) + 19) / 720,
                                  (NN**4 + 10*NN**2 + 60*NN*kn**3 - 60*NN*kn -
                                   30*kn**4 + kn**2*(60 - 30*NN**2) - 11) / 720,
                                  (NN**4 - 20*NN**2 - 30*kn**4 + kn**3*(60*NN + 120) + kn**2*(-30*NN**2 - 180*NN - 120) + kn*(60*NN**2 + 120*NN) + 19) / 720]
                    elif pp == 5:
                        T_p_mu = [(5*NN**4 - 250*NN**2 + 720*NN - 12*kn**5 + kn**4*(30*NN - 150) + kn**3*(-20*NN**2 + 300*NN - 700) + kn**2*(-150*NN**2 + 1050*NN - 1500) + kn*(2*NN**4 - 350*NN**2 + 1500*NN - 1440) - 475) / 1440,
                                  (3*NN**4 - 30*NN**2 - 12*kn**5 + kn**4*(30*NN - 90) + kn**3*(-20*NN**2 + 180*NN - 220) +
                                   kn**2*(-90*NN**2 + 330*NN - 180) + kn*(2*NN**4 - 110*NN**2 + 180*NN) + 27) / 1440,
                                  (NN**4 + 10*NN**2 - 12*kn**5 + kn**4*(30*NN - 30) + kn**3*(-20*NN**2 + 60*NN + 20) +
                                   kn**2*(-30*NN**2 - 30*NN + 60) + kn*(2*NN**4 + 10*NN**2 - 60*NN) - 11) / 1440,
                                  (-NN**4 - 10*NN**2 - 12*kn**5 + kn**4*(30*NN + 30) + kn**3*(-20*NN**2 - 60*NN + 20) + kn**2*(
                                      30*NN**2 - 30*NN - 60) + kn*(2*NN**4 + 10*NN**2 + 60*NN) + 11) / 1440,
                                  (-3*NN**4 + 30*NN**2 - 12*kn**5 + kn**4*(30*NN + 90) + kn**3*(-20*NN**2 - 180*NN - 220) + kn**2*(90*NN**2 + 330*NN + 180) + kn*(2*NN**4 - 110*NN**2 - 180*NN) - 27) / 1440]
                    else:
                        T_p_mu = make_polynom_table2(kn, pp, varN=NN)

                    num_coeff = sym.Poly(rem, qq).all_coeffs()[-1::-1]
                    if len(num_coeff) > len(T_p_mu):
                        raise ValueError(
                            "Error: polynomial division may be wrong!")
                    # Sum up
                    result = 0
                    for i in range(len(num_coeff)):
                        result += num_coeff[i] * T_p_mu[i] * \
                            sym.exp(j * 2 * pi * n0 / self.N * i * sg)
                    result *= sym.exp(-j * 2 * pi * n0 /
                                      self.N * k * sg) / c0**pp

                    return const * (result_0 + result)
        return None

    def term(self, expr, n, k, **kwargs):

        const, expr = factor_const(expr, n)
        if self.is_inverse:
            const /= self.N

        if expr.has(AppliedUndef):
            # Handle v(n), v(n) * y(n), 3 * v(n) / n etc.
            result = const * self.function(expr, n, k)
            return result

        # Transforms of type x[n] o--o X(q) with q=exp(-j*2*pi/N*k)
        # (mainly used with the forward transform, but works also for
        # backward transforms)

        # Call transform, use inversion for first case (change sign
        # later back for second case)
        k = -k if self.is_inverse else k

        q = sym.Symbol('q')
        res = self.termXq(expr, n, k, q, 0, self.N - 1)
        # Find handles
        if not res is None:
            # Simplify q**N terms
            res.simp_qN(q, self.N)
            # Make final transform by putting all cases and terms together
            res.make_transform(self.N, self.is_inverse, q,
                               sym.exp(-sym.I * 2 * pi / self.N * k), k,
                               kwargs.get('piecewise', None))
            result = res.Xk

            # TODO smart simplification of result since simplify takes
            # too long or makes sometimes weird expressions might be
            # better in make_transform

            # result = sym.simplify(result)   #only useful if handles found
            return const * result

        # change back, as next handles do not need  replacement k -> -k
        if self.is_inverse:
            k = -k

        # Transforms of type x[n] = rational function(qq)  with qq=exp(j*2*pi/N*n)
        # (mainly used with the backward transform, but works also for forward transform)
        result = self.termXk(expr, n, k)
        if not result is None:
            return const * result

        # No special handle identified, general expression
        k = -k if self.is_inverse else k
        if not self.N.is_symbol:
            result = sym.summation(
                expr * sym.exp(-j * 2 * pi * n * k / self.N), (n, 0, self.N - 1))
        else:
            result = sym.Sum(expr * sym.exp(-j * 2 * pi *
                             n * k / self.N), (n, 0, self.N - 1))
        return const * result


dft_transformer = DFTTransformer()


def discrete_fourier_transform(expr, n, k, N=None, evaluate=True,
                               **kwargs):
    """Compute bilateral discrete Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(k)
    """

    return dft_transformer.transform(expr, n, k, evaluate=evaluate, N=N,
                                     **kwargs)


def DFT(expr, n, k, N=None, evaluate=True, **kwargs):
    """Compute bilateral discrete Fourier transform of expr.

    Undefined functions such as x(n) are converted to X(k)
    """

    return dft_transformer.transform(expr, n, k, evaluate=evaluate, N=N,
                                     **kwargs)


def DFTmatrix(N):
    """Return DFT matrix of size `N` x `N`."""

    from .functions import exp
    from .sym import j, pi

    w = exp(-j * 2 * pi / N)

    a = Matrix.zeros(N)

    for row in range(N):
        for col in range(N):
            a[row, col] = w ** (row * col)
    return a
