"""
This module supports simple linear circuit analysis based on the
following ideal components:

V independent voltage source
I independent current source
R resistor
C capacitor
L inductor

These components are converted to s-domain models and so capacitor and
inductor components can be specified with initial voltage and
currents, respectively, to model transient responses.

The components are represented by either Thevenin or Norton one-port
networks with the following attributes:

Zoc open-circuit impedance
Ysc short-circuit admittance
Voc open-circuit voltage
Isc short-circuit current

Components can either be connected in series (+) or parallel (|).

Components can be connected to form two-port networks.  Methods are
provided to determine transfer responses.

Note, the currents or voltages may cancel and thus a Norton or
Thevenin network may collapse to a pure impedance.  However, the
impedance cannot cancel causing the collapse to a voltage or current
source.

To print the rational functions in canonical form (with the highest
power of s in the denominator with a unity coefficient), use
print(x.canonical()) or x.pprint() for pretty printing.

Some multiport networks, such as a shunt R, have a singular Z matrix.
Thus switching to the Y matrix and back to the Z matrix produces a
bogus result.  The same thing occurs for a series R; this has a
singular Y matrix.

Perhaps could have Ser and Par classes that store the underlying
one-ports so that the network can be reconstructed; similar to how
sympy stores expressions?  These could have internal lists of all the
networks that are in series or parallel.

It also may be better to defer the choice of the two-port model.  For
example, a T section would store the three sub-networks rather than
generating a B matrix.  The appropriate model would be generated when
desired.  This would avoid the inversion of singular matrices. The
downside is that each object would require methods to generate each
type of two-port model.

The original implementation stored _Expr as rational functions (using
MRF class).  This is much faster, avoids symbolic inverse Laplace
transforms, but does not handle delays and symbolic component values.

TODO: Perhaps don't store R, V, I, as one-ports?  Or perhaps have a
source-free one-port class?

TODO: Fix handling of buffered two ports (amplifier / delay).

Copyright 2014 Michael Hayes, UCECE
"""

# Consider chaining a resistor shunt to a two-port network described by the
# A matrix A1.  The shunt has an A matrix A2.  The result has an A matrix A3.
#
# A1 = [A11 A12]
#      [A21 A22]
#
# A2 = [1     0]
#      [1/R   1]
#
# A3 = A1 * A2 = [A11 + A12/R  A12]
#                [A21 + A22/R  A22]
#
# We should get the same result by adding the Y matrices.
#
# Y1 = [A22/A12  -A11 A22 / A12 + A21]
#      [-1/A12                A11/A12]
#
# Y3 = [A22/A12  -A11 A22 / A12 + A21]
#      [-1/A12          A11/A12 + 1/R]
#
# but
#
# Y2 = [inf   inf]
#      [inf   inf]
#
# We can get the correct answer using:
#
# A2 = lim x->0  [1    x]
#                [1/R  1]
#
# Now det(A2) = lim x->0 (1 - x/R)
#
# and Y2 = lim x->0 [1/x   (1/x - 1/R)]
#                   [-1/x          1/x]
#
# Note when x=0, then
#
# Y2 = lim x->0 [1/0   1/0]
#               [-1/0  1/0]       
#
# and we lose the information on R.
#
# The same problem occurs with a series R and the Z matrix. 
#
# A2 = lim x->0  [1    R]
#                [x    1]
#
# Now det(A2) = lim x->0 (1 - R x)
#
# and Z2 = lim x->0 [1/x  1/x - R]
#                   [1/x      1/x]
#
# Thus it is advantageous to represent two-ports by the A (or B)
# matrix.  However, things will go wrong when we transform to the Y or
# Z matrix for specific cases.

from __future__ import division
from warnings import warn
import numpy as np
import sympy as sym
from sympy.utilities.lambdify import lambdify


def _mrffmt(prefix, arg):
    
    prefix += ' = '
    str = arg.__str__()
    pad = ' ' * len(prefix)
    return prefix + str.replace('\n', '\n' + pad)
    

def _strpair(prefix1, arg1, prefix2, arg2):

    return '%s\n\n%s' % (_mrffmt(prefix1, arg1), _mrffmt(prefix2, arg2))


def _pretty_fmt(prefix, arg):
    
    if hasattr(arg, 'val'):
        arg = arg.val

    return sym.pretty(sym.Eq(sym.Symbol(prefix), arg))


def _pretty_strpair(prefix1, arg1, prefix2, arg2):

    return '%s\n\n%s' % (_pretty_fmt(prefix1, arg1), _pretty_fmt(prefix2, arg2))


def pprint(expr):

    print(pretty(expr))


def pretty(expr):

    if hasattr(expr, 'pretty'):
        return expr.pretty()
    else:
        return sym.pretty(expr)


def DeltaWye(Z1, Z2, Z3):

    ZZ = (Z1 * Z2 + Z2 * Z3 + Z3 * Z1)
    return (ZZ / Z1, ZZ / Z2, ZZ / Z3)


def WyeDelta(Z1, Z2, Z3):

    ZZ = Z1 + Z2 + Z3
    return (Z2 * Z3 / ZZ, Z1 * Z3 / ZZ, Z1 * Z2 / ZZ)


def poles(expr, var=None):

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    numer, denom = expr.as_numer_denom()
    poles = sym.roots(sym.Poly(denom, var))
    return poles


def zeros(expr, var=None):

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    numer, denom = expr.as_numer_denom()
    zeros = sym.roots(sym.Poly(numer, var))
    return zeros


def residue(expr, var, pole, poles):

    # Remove pole from list of poles; sym.cancel
    # doesn't always work, for example, for complex poles.
    poles2 = poles.copy()
    poles2[pole] -= 1

    numer, denom = expr.as_numer_denom()
    D = sym.Poly(denom, var)
    K = D.LC()

    D = [(var - p) ** poles2[p] for p in poles2]
    denom = sym.Mul(K, *D)

    d = sym.limit(denom, var, pole)

    if d != 0:
        tmp = numer / denom
        return sym.limit(tmp, var, pole)

    print("Trying l'hopital's rule")
    tmp = numer / denom
    tmp = sym.diff(tmp, var)
    
    return sym.limit(tmp, var, pole)


def residues(expr, var=None):

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    N, D = _as_ratfun_parts(expr, var)

    Q, M = N.div(D)

    expr = M / D

    P = poles(expr, var)
    F = []
    R = []
    for p in P:
        
        # Number of occurrences of the pole.
        N = P[p]

        f = var - p

        if N == 1:
            F.append(f)
            R.append(residue(expr, var, p, P))
            continue

        # Handle repeated poles.
        expr2 = expr * f ** N
        for n in range(1, N + 1):
            m = N - n
            F.append(f ** n)
            R.append(sym.limit(sym.diff(expr2, var, m), var, p) / sym.factorial(m))

    return F, R, Q


def partfrac(expr, var=None):
    """Convert rational function into partial fraction form (standard form).

    See also canonical, general, and ZPK"""

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    ratfun, delay = _as_ratfun_delay(expr, var)

    F, R, Q = residues(ratfun, var)

    expr = Q
    for f, r in zip(F, R):
        expr = expr + r / f

    if delay != 0:
        expr *= sym.exp(-var * delay)

    return expr


def _as_ratfun_parts(expr, var=None):
    
    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    if not expr.is_rational_function():
        raise ValueError('Expression not a rational function')

    numer, denom = expr.as_numer_denom()
    N = sym.Poly(numer, var)
    D = sym.Poly(denom, var)
    
    return N, D


def _as_ratfun_delay(expr, var=None):
    
    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    F = expr.as_ordered_factors()

    delay = sym.sympify(0)
    ratfun = sym.sympify(1)
    for f in F:
        b, e = f.as_base_exp()
        if b == sym.E and e.is_polynomial(var):
            p = sym.Poly(e, var)
            c = p.all_coeffs()
            if p.degree() == 1:
                delay -= c[0]
                if c[1] != 0:
                    ratfun *= sym.exp(c[1])
                continue

        ratfun *= f
            
    if not ratfun.is_rational_function():
        raise ValueError('Expression not a product of rational function and exponential')

    return ratfun, delay


def general(expr, var=None):
    """Convert rational function into general form.

    See also canonical, partfrac, and ZPK"""

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    return sym.cancel(expr, var)



def canonical(expr, var=None):
    """Convert rational function into canonical form.

    See also general, partfrac, and ZPK"""


    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    ratfun, delay = _as_ratfun_delay(expr, var)

    N, D = _as_ratfun_parts(ratfun, var)

    K = sym.cancel(N.LC() / D.LC())
    sym.pprint(K)
    if delay != 0:
        K = K * sym.exp(var * delay)

    # Divide by leading coefficient
    N = N.monic()
    D = D.monic()
    
    return K * (N / D)


def ZPK(expr, var=None):
    """Convert rational function into ZPK form.

    See also canonical, general, and partfrac"""

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    ratfun, delay = _as_ratfun_delay(expr, var)

    N, D = _as_ratfun_parts(ratfun, var)

    K = sym.cancel(N.LC() / D.LC())
    if delay != 0:
        K = K * sym.exp(var * delay)

    zeros = sym.roots(N)
    zz = [(var - z) for z in zeros]

    poles = sym.roots(D)
    pp = [1 / (var - p) for p in poles]
        
    return sym.Mul(K, *(zz + pp))


def _inverse_laplace(expr, var, t):

    ratfun, delay = _as_ratfun_delay(expr, var)

    N, D = _as_ratfun_parts(ratfun, var)

    Q, M = N.div(D)

    result1 = 0

    # Delayed time.
    td = t - delay

    if Q:
        print('Warning: Impulse response has impulses and/or derivatives of impulses.')
        C = Q.all_coeffs()
        for n, c in enumerate(C):
            result1 += c * sym.diff(sym.DiracDelta(td), t, len(C) - n - 1)
        
    expr = M / D

    P = poles(expr, var)
    result2 = 0
    for p in P:
        
        # Number of occurrences of the pole.
        N = P[p]

        f = var - p

        if N == 1:
            r = residue(expr, var, p, P)
            result2 += r * sym.exp(p * td)
            continue

        # Handle repeated poles.
        expr2 = expr * f ** N
        for n in range(1, N + 1):
            m = N - n
            r = sym.limit(sym.diff(expr2, var, m), var, p) / sym.factorial(m)
            result2 += r * sym.exp(p * td) * td**(n - 1)

    return result1 + result2 * sym.Heaviside(td)


def inverse_laplace(expr, s, t):
    """Determine inverse Laplace transform of expression"""

    try:
        result = _inverse_laplace(expr, s, t)
    except:
        
        # Try splitting into partial fractions to help sympy.
        expr = partfrac(expr, s)

        # This barfs when needing to generate Dirac deltas
        from sympy.integrals.transforms import inverse_laplace_transform
        result = inverse_laplace_transform(expr, s, t)

    return result


def initial_value(expr, var=None):

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    return sym.limit(expr * var, var, sym.oo)


def final_value(expr, var=None):

    if hasattr(expr, 'expr'):
        var = expr.s
        expr = expr.expr

    return sym.limit(expr * var, var, 0)


class _Expr(object):
    
    s, t, f = sym.symbols('s t f')
    
    
    @property
    def expr(self):    
        return self.val
    
    
    def __init__(self, val, simplify=True, real=False):
        
        if isinstance(val, _Expr):
            val = val.val
            
        if real and isinstance(val, str):
            val = sym.symbols(val, real=True)

        if simplify:
            val = sym.sympify(val).cancel()

        self.val = val


    def __str__(self):

        return self.val.__str__()


    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.val.__str__())

    
    def __neg__(self):
        """Negation"""
        
        return self.__class__(-self.val)


    def __rdiv__(self, x):
        """Reverse divide"""

        x = _Expr(x)
        return self.__class__(x.val / self.val)


    def __rtruediv__(self, x):
        """Reverse true divide"""
            
        x = _Expr(x)
        return self.__class__(x.val / self.val)


    def __mul__(self, x):
        """Multiply"""
        
        x = _Expr(x)
        return self.__class__(self.val * x.val)


    def __rmul__(self, x):
        """Reverse multiply"""
        
        x = _Expr(x)
        return self.__class__(self.val * x.val)


    def __div__(self, x):
        """Divide"""

        x = _Expr(x)
        return self.__class__(self.val / x.val)


    def __truediv__(self, x):
        """True divide"""

        x = _Expr(x)
        return self.__class__(self.val / x.val)
    

    def __add__(self, x):
        """Add"""
        
        x = _Expr(x)
        return self.__class__(self.val + x.val)
    

    def __radd__(self, x):
        """Reverse add"""
        
        x = _Expr(x)
        return self.__class__(self.val + x.val)
    
    
    def __rsub__(self, x):
        """Reverse subtract"""
        
        x = _Expr(x)
        return self.__class__(x.val - self.val)
    

    def __sub__(self, x):
        """Subtract"""
        
        x = _Expr(x)
        return self.__class__(self.val - x.val)
    
    
    def __or__(self, x):
        """Parallel combination"""
        
        return self.parallel(x)
    
    
    def __eq__(self, x):
        """Equality"""

        x = _Expr(x)
        return self.val == x.val


    def __ne__(self, x):
        """Inequality"""

        x = _Expr(x)
        return self.val != x.val


    def parallel(self, x):
        """Parallel combination"""
        
        x = _Expr(x)
        return self.__class__(self.val * x.val / (self.val + x.val))
    
    
    def differentiate(self):
        """Differentiate (multiply by s)"""
        
        return self.__class__(self.val * self.s)
    

    def integrate(self):
        """Integrate (divide by s)"""
        
        return self.__class__(self.val / self.s)
    
    
    def delay(self, T):
        """Apply delay of T seconds by multiplying by exp(-s T)"""
        
        T = _Expr(T)
        return self.__class__(self.val * sym.exp(-T * self.s))


    def zeros(self):
        
        return zeros(self.expr, self.s)


    def poles(self):
        
        return poles(self.expr, self.s)


    def residues(self):
        
        return residues(self.expr, self.s)


    @property
    def canonical(self):
        """Convert rational function to canonical form with unity
        highest power of denominator.

        See also general, partfrac, and ZPK"""
        
        return self.__class__(canonical(self.expr, self.s), simplify=False)


    @property
    def general(self):
        """Convert rational function to general form

        See also canonical, partfrac, and ZPK"""
        
        return self.__class__(general(self.expr, self.s), simplify=False)


    @property
    def partfrac(self):
        """Convert rational function into partial fraction form.

        See also canonical, general, and ZPK"""

        return self.__class__(partfrac(self.expr, self.s), simplify=False)


    @property
    def ZPK(self):
        """Convert to pole-zero-gain (PZK) form.
        
        See also canonical, general, and partfrac"""
        
        return self.__class__(ZPK(self.expr, self.s), simplify=False)


    @property
    def initial_value(self):
        """Determine value at t = 0"""
        
        return initial_value(self.expr, self.s)


    @property
    def final_value(self):
        """Determine value at t = oo"""
        
        return final_value(self.expr, self.s)


    def _as_ratfun_parts(self):
        
        return _as_ratfun_parts(self.expr, self.s)


    def split_strictly_proper(self):
        
        N, D = self._as_ratfun_parts()

        Q, M = N.div(D)

        return Q.as_expr(), M / D


    def mrf(self):

        from mrf import MRF

        val = self.val.cancel()

        N, D = self._as_ratfun_parts()

        from scipy import poly1d
        Np = poly1d([float(x) for x in N.all_coeffs()])
        Dp = poly1d([float(x) for x in D.all_coeffs()])
        
        return MRF(Np.c, Dp.c)


    def pprint(self):
        """Pretty print"""
        sym.pprint(self.val)


    def pretty(self):
        """Make pretty string"""
        return sym.pretty(self.val)


    def simplify(self):
        """Simplify"""
        # Not sure if this is needed; perhaps use canonical?
        val = sym.ratsim(self.val)
        return self.__class__(val)
    
    
    def inverse_laplace(self):
        """Attempt inverse Laplace transform"""
        
        print('Determining inverse Laplace transform...')
        return inverse_laplace(self.expr, self.s, self.t)


    def transient_response(self, t=None):
        """Evaluate transient (impulse) response"""
        
        expr = self.inverse_laplace()
        if t == None:
            return expr

        print('Evaluating inverse Laplace transform...')
        
        func = lambdify(self.t, expr, ("numpy", "sympy", "math"))
        
        try:
            # FIXME for scalar t
            response = np.array([complex(func(t1)) for t1 in t])
            # The following does not work if all the sympy functions are not
            # converted to numpy functions.
            # response = func(t)
            
        except NameError:
            raise RuntimeError('Cannot compute inverse Laplace transform')
        
        # The result should be real so quietly remove any imaginary
        # component.
        return response.real
    
    
    def impulse_response(self, t=None):
        """Evaluate transient (impulse) response"""
        
        return self.transient_response(t)


    def step_response(self, t=None):
        """Evaluate step response"""
        
        return (self / self.s).transient_response(t)
    
    
    def frequency_response(self, f=None):
        """Evaluate frequency response"""
        
        expr = self.val.subs(self.s, sym.I * 2 * sym.pi * self.f)
        
        if f == None:
            return expr
        
        func = lambdify(self.f, expr, modules="numpy")
        return np.array([func(f1) for f1 in f])



class Zs(_Expr):
    """s-domain impedance value"""

    @classmethod
    def C(cls, Cval):
    
        Cval = _Expr(Cval, real=True)
        return cls(1 / Cval).integrate()


    @classmethod
    def G(cls, Gval):
    
        Gval = _Expr(Gval, real=True)
        return cls(1 / Gval)


    @classmethod
    def L(cls, Lval):
    
        Lval = _Expr(Lval, real=True)
        return cls(Lval).differentiate()


    @classmethod
    def R(cls, Rval):
    
        Rval = _Expr(Rval, real=True)
        return cls(Rval)


class Ys(_Expr):
    """s-domain admittance value"""
    

    @classmethod
    def C(cls, Cval):
    
        Cval = _Expr(Cval, real=True)
        return cls(Cval).differentiate()


    @classmethod
    def G(cls, Gval):
    
        Gval = _Expr(Gval, real=True)
        return cls(Gval)


    @classmethod
    def L(cls, Lval):
    
        Lval = _Expr(Lval, real=True)
        return cls(1 / Lval).integrate()


    @classmethod
    def R(cls, Rval):
    
        Rval = _Expr(Rval, real=True)
        return cls(1 / Rval)


class Vs(_Expr):
    """s-domain voltage (units V s / radian)"""
    pass


class Is(_Expr):
    """s-domain current (units A s / radian)"""
    pass


class Avs(_Expr):
    """s-domain voltage ratio"""
    pass


class Ais(_Expr):
    """s-domain current ratio"""
    pass


class OnePort(object):

    # Attributes: Y, Z, V, I

    def __add__(self, x):
        """Series combination"""

        return self.series(x)


    def __or__(self, x):
        """Parallel combination"""

        return self.parallel(x)


    @property
    def Zoc(self):    
        return self.Z


    @property
    def Voc(self):    
        return self.V


    @property
    def Ysc(self):    
        return self.Y


    @property
    def Isc(self):    
        return self.I


    def ladder(self, *args):
        """Create (unbalanced) ladder network"""

        return Ladder(self, *args)


    def lsection(self, OP2):
        """Create L section (voltage divider)"""

        if not issubclass(OP2.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        return LSection(self, OP2)


    def tsection(self, OP2, OP3):
        """Create T section"""

        if not issubclass(OP2.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        if not issubclass(OP3.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        return TSection(self, OP2, OP3)


class Norton(OnePort):
    """Norton (Y) model

            +-------------------+  
    I1      |                   |      -I1
    -->-+---+        Y          +---+--<--
        |   |                   |   |
        |   +-------------------+   |
        |                           |
        |                           |
        |          +------+         |
        |          |      |         |
        +----------+ -I-> +---------+
                   |      |    
                   +------+    

     +                V1                 -
    """

    def __init__(self, Yval, Ival=Is(0)):

        #print('<N> Y:', Yval, 'I:', Ival)
        if not isinstance(Yval, Ys):
            raise ValueError('Yval not Ys')
        if not isinstance(Ival, Is):
            raise ValueError('Ival not Is')
        self.Y = Yval
        self.I = Ival


    @property
    def Z(self):    
        return Zs(1 / self.Y)


    @property
    def V(self):    
        return Vs(self.I / self.Y)


    def __repr__(self):

        return '%s(%s, %s)' % (self.__class__.__name__, self.Y.__repr__(), self.I.__repr__())


    def __str__(self):
        
        return _strpair('Y', self.Y, 'I', self.I)


    def pretty(self):

        return _pretty_strpair('Y', self.Y, 'I', self.I)


    def thevenin(self):

        return Thevenin(self.Z, self.V)


    def norton(self):

        return self


    def series(self, x):

        return self.thevenin().series(x).norton()


    def parallel(self, x):

        if isinstance(x, (C, G, L, R, V, Thevenin)):
            if x.Z == 0:
                raise ValueError('Cannot connect voltage source in parallel.')
            y = x.norton()
            return Norton(self.Y + y.Y, self.I + y.I)
        elif isinstance(x, (I, Norton)):
            return Norton(self.Y + x.Y, self.I + x.I)
        else:
            raise ValueError('Unhandled type ', type(x))            


class Thevenin(OnePort):
    """Thevenin (Z) model

         +------+    +-------------------+  
    I1   | +  - |    |                   | -I1
    -->--+  V   +----+        Z          +--<--
         |      |    |                   |    
         +------+    +-------------------+  

    +                       V1                -
    """


    def __init__(self, Zval, Vval=Vs(0)):

        #print('<T> Z:', Zval, 'V:', Vval)
        if not isinstance(Zval, Zs):
            raise ValueError('Zval not Zs')
        if not isinstance(Vval, Vs):
            raise ValueError('Vval not Vs')
        self.Z = Zval
        self.V = Vval


    @property
    def Y(self):    
        return Ys(1 / self.Z)


    @property
    def I(self):    
        return Is(self.V / self.Z)


    def __repr__(self):

        return '%s(%s, %s)' % (self.__class__.__name__, self.Z.__repr__(), self.V.__repr__())


    def __str__(self):

        return _strpair('Z', self.Z, 'V', self.V)


    def pretty(self):

        return _pretty_strpair('Y', self.Y, 'I', self.I)


    def norton(self):

        return Norton(self.Y, Is(self.V / self.Z))


    def thevenin(self):

        return self


    def series(self, x):

        if isinstance(x, (C, G, L, R, V, Thevenin)):
            return Thevenin(self.Z + x.Z, self.V + x.V)
        elif isinstance(x, (I, Norton)):
            if x.Y == 0:
                raise ValueError('Cannot connect current source in series.')
            y = x.thevenin()
            return Thevenin(self.Z + y.Z, self.V + y.V)
        else:
            raise ValueError('Unhandled type ', type(x))            


    def parallel(self, x):

        return self.norton().parallel(x).thevenin()


    def parallel_ladder(self, *args):
        """Add unbalanced ladder network in parallel; alternately in parallel and series.

          +---------+       +---------+       
       +--+   self  +---+---+   Z1    +---+---
       |  +---------+   |   +---------+   |
       |              +-+-+             +-+-+
       |              |   |             |   |
       |              |Z0 |             |Z2 |
       |              |   |             |   |
       |              +-+-+             +-+-+                       
       |                |                 |
       +----------------+-----------------+---
       """

        ret = self
        for m, arg in enumerate(args):
            if m & 1:
                ret = ret.series(arg)
            else:
                ret = ret.parallel(arg)
        return ret


    def parallel_C(self, Z0, Z1, Z2):
        """Add C network in parallel.

          +---------+      +---------+        
       +--+   self  +------+   Z0    +---+----
       |  +---------+      +---------+   |    
       |                               +-+-+
       |                               |   |
       |                               |Z1 |
       |                               |   |
       |                               +-+-+                       
       |                   +---------+   |
       +-------------------+   Z2    +---+----
                           +---------+
       """

        return self.series(Z0).series(Z2).parallel(Z1)



    def parallel_L(self, Z0, Z1):
        """Add L network in parallel.

          +---------+      +---------+        
       +--+   self  +------+   Z0    +---+----
       |  +---------+      +---------+   |    
       |                               +-+-+
       |                               |   |
       |                               |Z1 |
       |                               |   |
       |                               +-+-+                       
       |                                 |
       +---------------------------------+----
       """

        return self.series(Z0).parallel(Z1)

                   
    def parallel_pi(self, Z0, Z1, Z2):
        """Add Pi (Delta) network in parallel.

          +---------+       +---------+       
       +--+   self  +---+---+   Z1    +---+---
       |  +---------+   |   +---------+   |
       |              +-+-+             +-+-+
       |              |   |             |   |
       |              |Z0 |             |Z2 |
       |              |   |             |   |
       |              +-+-+             +-+-+                       
       |                |                 |
       +----------------+-----------------+---
       """
                   
        return (self.parallel(Z0) + Z1).parallel(Z2)


    def parallel_T(self, Z0, Z1, Z2):
        """Add T (Y) network in parallel.

          +---------+       +---------+        +---------+       
       +--+   self  +-------+   Z0    +---+----+   Z2    +---
       |  +---------+       +---------+   |    +---------+   
       |                                +-+-+
       |                                |   |
       |                                |Z1 |
       |                                |   |
       |                                +-+-+                       
       |                                  |
       +----------------------------------+------------------
       """

        return (self.parallel(Z0) + Z1).parallel(Z2)


    def load(self, x):
        """Apply a load and create a Load object that stores the voltage
        across the load and the current through it"""
        
        # This may need some pondering.  What if a Thevenin network is
        # connected?
        return Load(self.parallel(x).V, self.series(x).I)


class Load(object):

    def __init__(self, Vval, Ival):

        self.V = Vval
        self.I = Ival


class R(Thevenin):
    """Resistor"""

    @property
    def Z(self):
        return Zs.R(self.R)


    @property
    def V(self):
        return Vs(0)


    def __init__(self, Rval):
    
        self.R = _Expr(Rval)


    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.R.__str__())


class G(Norton):
    """Conductance"""

    @property
    def Z(self):
        return Zs.G(self.G)


    @property
    def V(self):
        return Vs(0)


    def __init__(self, Gval):
    
        self.G = _Expr(Gval)


    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.G.__str__())


class L(Thevenin):
    """Inductor"""


    @property
    def Z(self):
        return Zs.L(self.L)

    
    @property
    def V(self):
        return Vs(-self.i0 * self.L)


    def __init__(self, Lval, i0=0.0):
    
        self.L = _Expr(Lval)
        self.i0 = _Expr(i0)


    def __repr__(self):

        if self.i0 != 0:
            return '%s(%s, %s)' % (self.__class__.__name__, self.L.__str__(), self.i0.__str__())
        else:
            return '%s(%s)' % (self.__class__.__name__, self.L.__str__())


class C(Thevenin):
    """Capacitor"""


    @property
    def Z(self):
        return Zs.C(self.C)

    
    @property
    def V(self):
        return Vs(self.v0).integrate()
        

    def __init__(self, Cval, v0=0.0):
    
        self.C = _Expr(Cval)
        self.v0 = _Expr(v0)


    def __repr__(self):
        
        if self.v0 != 0:
            return '%s(%s, %s)' % (self.__class__.__name__, self.C.__str__(), self.v0.__str__())
        else:
            return '%s(%s)' % (self.__class__.__name__, self.C.__str__())



class Xtal(Thevenin):
    """Crystal"""

    def __init__(self, C0, R1, L1, C1):
        """Create a crystal using a Butterworth van Dyke model;
        harmonic resonances are not modelled

        This is modelled as a series R, L, C circuit in parallel
        with C0."""
    
        xtal = (R(R1) + L(L1) + C(C1)) | C(C0)

        super (Xtal, self).__init__(xtal.Z, xtal.V)


class FerriteBead(Thevenin):
    """Ferrite bead"""

    def __init__(self, Rs, Rp, Cp, Lp):
        """Create a ferrite bead (lossy inductor)

        This is modelled as a series resistor (Rs) connected 
        to a parallel R, L, C network (Rp, Lp, Cp).
        """

        ferrite = R(Rs) + (R(Rp) + L(Lp) + C(Cp))

        super (Xtal, self).__init__(ferrite.Z, xtal.V)


class V(Thevenin):
    """Voltage source (note a voltage source of voltage V1 has
    a s domain voltage of V1 / s."""

    @property
    def Z(self):
        return Zs(0)

    
    @property
    def V(self):
        return Vs(self.v.integrate())


    def __init__(self, v):
        """Create voltage source with voltage v V."""
    
        self.v = _Expr(v)


    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.v.__str__())


    def __str__(self):
        
        return self.V.__str__()


    def series(self, x):

        if isinstance(x, V):
            return V(self.V + x.V)
        return super (V, self).series(x)


    def parallel(self, x):

        raise ValueError('Cannot connect voltage source in parallel.')


    def thevenin(self):

        return Thevenin(Zs(0), self.V)


    def norton(self):

        warn('Converting a voltage source to a Norton network is dodgy...')
        return Norton(Ys(1 / Zs(0)), self.I)


class I(Norton):
    """Current source (note a current source of current I1 has
    a s domain current of I1 / s."""

    @property
    def Y(self):
        return Ys(0)

    
    @property
    def I(self):
        return Is(self.i.integrate())


    def __init__(self, i):
        """Create current source with current i A."""
    
        self.i = _Expr(i)


    def __repr__(self):

        return '%s(%s)' % (self.__class__.__name__, self.i.__str__())


    def __str__(self):
        
        return self.I.__str__()


    def parallel(self, x):

        if isinstance(x, I):
            return I(self.I + x.I)
        return super (I, self).parallel(x)


    def series(self, x):

        raise ValueError('Cannot connect current source in series.')


    def thevenin(self):

        warn('Converting a current source to a Thevenin network is dodgy...')
        return Thevenin(Zs(1 / Ys(0)), self.V)


    def norton(self):

        return Norton(Ys(0), self.I)



class Vsector(sym.Matrix):

    # Unlike numpy.ndarray, the sympy.Matrix runs all the elements
    # through sympify, creating sympy objects and thus losing the
    # original type information and associated methods.

    def __new__ (cls, *args):

        if len(args) == 2:
            return super (Vsector, cls).__new__(cls, (args[0], args[1]))

        return super (Vsector, cls).__new__(cls, *args)


    def __getitem__(self, key):

        item = super (Vsector, self).__getitem__(key)

        if isinstance(key, int):
            return self._typewrap(item)
        warn('Returning sympy object')
        return item


class VsVector(Vsector):
    
    _typewrap = Vs


class IsVector(Vsector):
    
    _typewrap = Is


class YVector(Vsector):

    _typewrap = Ys


class ZVector(Vsector):

    _typewrap = Zs


class _TwoPortMatrix(sym.Matrix):

    def __new__ (cls, *args):

        if len(args) == 4:
            return super (_TwoPortMatrix, cls).__new__(cls, ((args[0], args[1]), (args[2], args[3])))

        return super (_TwoPortMatrix, cls).__new__(cls, *args)

    
    # The following properties are fallbacks when other conversions have
    # not been defined.

    @property
    def A(self):
        return AMatrix(self.B.inv())


    @property
    def B(self):
        return BMatrix(self.A.inv())


    @property
    def G(self):
        return GMatrix(self.H.inv())


    @property
    def H(self):
        return HMatrix(self.G.inv())


    @property
    def Y(self):
        return YMatrix(self.Z.inv())


    @property
    def Z(self):
        return ZMatrix(self.Y.inv())



class AMatrix(_TwoPortMatrix):
    """
    +-  -+     +-       -+   +-  -+    
    | V1 |  =  | A11  A12|   | V2 | 
    | I1 |     | A21  A22|   |-I2 | 
    +-  -+     +-       -+   +-  -+    

           +-         -+
    units  | 1     ohm |
           | 1/ohm   1 |
           +-         -+ 

    A buffered two-port has A12 = A22 = 0.

    A = inv(B)
    """

    @property
    def A11(self):
        return self[0, 0]


    @property
    def A12(self):
        return self[0, 1]


    @property
    def A21(self):
        return self[1, 0]


    @property
    def A22(self):
        return self[1, 1]


    @property
    def A(self):
        # Perhaps we should make a copy?
        return self


    @property
    def B(self):

        # Inverse
        det = self.det()
        if det == 0:
            warn('Producing dodgy B matrix')
        return BMatrix(self.A22 / det, -self.A12 / det, -self.A21 / det, self.A11 / det)


    @property
    def H(self):

        if self.A22 == 0:
            warn('Producing dodgy H matrix')
        return HMatrix(self.A12 / self.A22, self.det() / self.A22, -1 / self.A22, self.A21 / self.A22)


    @property
    def Y(self):

        # This produces a bogus Y matrix when A12 is zero (say for a
        # shunt element).   Note, it doesn't use A21.
        if self.A12 == 0:
            warn('Producing dodgy Y matrix')
        return YMatrix(self.A22 / self.A12, -self.det() / self.A12, -1 / self.A12, self.A11 / self.A12)


    @property
    def Z(self):

        # This produces a bogus Z matrix when A21 is zero (say for a
        # series element).   Note, it doesn't use A12.
        if self.A21 == 0:
            warn('Producing dodgy Z matrix')
        return ZMatrix(self.A11 / self.A21, self.det() / self.A21, 1 / self.A21, self.A22 / self.A21)


    @property
    def Z1oc(self):
        """open-circuit input impedance"""
        # Z11
        return Zs(self.A11 / self.A21)


    @classmethod
    def Zseries(cls, Zval):

        if not isinstance(Zval, Zs):
            raise ValueError('Zval not Zs')            

        return cls(1, Zval, 0, 1)


    @classmethod
    def Yseries(cls, Yval):

        if not isinstance(Yval, Ys):
            raise ValueError('Yval not Ys')

        return cls(1, 1 / Yval, 0, 1)


    @classmethod
    def Yshunt(cls, Yval):

        if not isinstance(Yval, Ys):
            raise ValueError('Yval not Ys')            

        return cls(1, 0, Yval, 1)


    @classmethod
    def Zshunt(cls, Zval):

        if not isinstance(Zval, Zs):
            raise ValueError('Zval not Zs')            

        return cls(1, 0, 1 / Zval, 1)


    @classmethod
    def transformer(cls, alpha):

        alpha = _Expr(alpha)

        return cls(1 / alpha, 0, 0, alpha)        


    @classmethod
    def gyrator(cls, R):

        R = _Expr(R)

        return cls(0, R, 1 / R, 0)        


    @classmethod
    def Lsection(cls, Z1, Z2):

        return cls.Zseries(Z1).chain(cls.Zshunt(Z2))


    @classmethod
    def Tsection(cls, Z1, Z2, Z3):

        return cls.Lsection(Z1, Z2).chain(cls.Zseries(Z3))


    @classmethod
    def Pisection(cls, Z1, Z2, Z3):

        return cls.Zshunt(Z1).chain(cls.Lsection(Z2, Z3))


    def chain(self, x):

        return self * x


    def cascade(self, x):

        return self.chain(x)



class BMatrix(_TwoPortMatrix):
    """
    +-  -+     +-       -+   +-  -+
    | V2 |  =  | B11  B12|   | V1 |
    |-I2 |     | B21  B22|   | I1 |
    +-  -+     +-       -+   +-  -+

           +-         -+
    units  | 1     ohm |
           | 1/ohm   1 |
           +-         -+ 


    B = inv(A)
    """

    @property
    def B11(self):
        return self[0, 0]


    @property
    def B12(self):
        return self[0, 1]


    @property
    def B21(self):
        return self[1, 0]


    @property
    def B22(self):
        return self[1, 1]


    @property
    def A(self):
        # Inverse
        det = self.det()
        return AMatrix(self.B22 / det, -self.B12 / det, -self.B21 / det, self.B11 / det)


    @property
    def B(self):
        # Perhaps we should make a copy?
        return self


    @property
    def G(self):

        return HMatrix(-self.B21 / self.B22, -1 / self.B22, self.det() / self.B22, -self.B12 / self.B22)


    @property
    def H(self):

        return HMatrix(-self.B12 / self.B11, 1 / self.B11, -self.det() / self.B11, -self.B21 / self.B11)


    @property
    def Y(self):

        return YMatrix(-self.B11 / self.B12, 1 / self.B12, self.det() / self.B12, -self.B22 / self.B12)


    @property
    def Z(self):

        return ZMatrix(-self.B22 / self.B21, -1 / self.B21, -self.det() / self.B21, -self.B11 / self.B21)


    @property
    def Z1oc(self):
        """open-circuit input impedance"""
        # Z11
        return Zs(-self.B22 / self.B21)


    @classmethod
    def Zseries(cls, Zval):

        if not isinstance(Zval, Zs):
            raise ValueError('Zval not Zs')            

        return cls(1, -Zval, 0, 1)


    @classmethod
    def Yseries(cls, Yval):

        if not isinstance(Yval, Ys):
            raise ValueError('Yval not Ys')

        return cls(1, -1 / Yval, 0, 1)


    @classmethod
    def Yshunt(cls, Yval):

        if not isinstance(Yval, Ys):
            raise ValueError('Yval not Ys')

        return cls(1, 0, -Yval, 1)


    @classmethod
    def Zshunt(cls, Zval):

        if not isinstance(Zval, Zs):
            raise ValueError('Zval not Zs')            

        return cls(1, 0, -1 / Zval, 1)


    @classmethod
    def voltage_amplifier(cls, Af, Ar=1e-9, Yin=1e-9, Zout=1e-9):
        """Voltage amplifier
        Af forward voltage gain
        Ar reverse voltage gain (ideally 0)
        Yin input admittance (ideally 0)
        Zout output impedance (ideally 0)
        """

        if Ar == 0 and Yin == 0 and Zout == 0:
            warn('Should use G matrix; tweaking B matrix to make invertible')
            Ar = 1e-9
            Yin = 1e-9
            Zout = 1e-9

        Af = _Expr(Af)
        Ar = _Expr(Ar)
        Yin = _Expr(Yin)
        Zout = _Expr(Zout)

        # This should be defined with a G matrix
        # 
        # G = [0   0]
        #     [Af  0]
        #
        # With this model, the inverse voltage gain is 1 / Af
        #
        # G = lim x->0  [0   x]
        #               [Af  0]
        # 
        # B = lim x->0  [Af    0/x]
        #               [0/x  -1/x]
        #
        # A = lim x->0  [1/Af 0/Af]
        #               [0/Af   -x]

        # Perhaps default Ar, Yin, and Zout to 1e-10 to get a reasonable
        # B matrix?

        return cls(1 / Ar, -1 / (Ar * Yin), -1 / (Ar * Zout), -1 / (Ar * Yin * Zout * (Af * Ar - 1)))


    @classmethod
    def current_amplifier(cls, Af, Ar=1e-9, Zin=1e-9, Yout=1e-9):
        """Current amplifier
        Af forward current gain
        Ar reverse current gain (ideally 0)
        Yin input admittance (ideally 0)
        Yout output impedance (ideally 0)
        """

        if Ar == 0 and Zin == 0 and Yout == 0:
            warn('Should use G matrix; tweaking B matrix to make invertible')
            Ar = 1e-9
            Zin = 1e-9
            Yout = 1e-9

        Af = _Expr(Af)
        Ar = _Expr(Ar)
        Zin = _Expr(Zin)
        Yout = _Expr(Yout)

        # This should be defined with a H matrix
        # 
        # H = [0   0]
        #     [Af  0]
        # 
        # With this model, the inverse current gain is 1 / Af
        #
        # H = lim x->0  [0   x]
        #               [Af  0]
        #
        # B = lim x->0  [1/x    0]
        #               [0    -Af]
        #
        # A = lim x->0  [1/x  -0/x]
        #               [0/x   -Af]

        return cls(1 / Ar, -1 / (Ar * Yout), -1 / (Ar * Zi), -1 / (Ar * Yo * Zi * (Af * Ar - 1)))


    @classmethod
    def voltage_differentiator(cls, Av=1):
        
        return cls.voltage_amplifier(_Expr(Av).differentiate())


    @classmethod
    def voltage_integrator(cls, Av):
        
        return cls.voltage_amplifier(_Expr(Av).integrate())


    @classmethod
    def current_differentiator(cls, Av):
        
        return cls.current_amplifier(_Expr(Av).differentiate())


    @classmethod
    def current_integrator(cls, Av):
        
        return cls.current_amplifier(_Expr(Av).integrate())


    @classmethod
    def transformer(cls, alpha):

        alpha = _Expr(alpha)

        return cls(alpha, 0, 0, 1 / alpha)        


    @classmethod
    def gyrator(cls, R):

        R = _Expr(R)

        return cls(0, R, 1 / R, 0)        


    @classmethod
    def Lsection(cls, Z1, Z2):

        Y = 1 / Z2
        return cls(1 + Y * Z1, -Z1, -Y, 1)        
        #return cls.Zseries(Z1).chain(cls.Zshunt(Z2))


    @classmethod
    def Tsection(cls, Z1, Z2, Z3):

        Y = 1 / Z2
        return cls(1 + Y * Z1, -Z1 -Z3 * (1 + Y * Z1), -Y, 1 + Y * Z3)        
        #return cls.Lsection(Z1, Z2).chain(cls.Zseries(Z3))


    @classmethod
    def Pisection(cls, Z1, Z2, Z3):

        return cls.Zshunt(Z1).chain(cls.Lsection(Z2, Z3))


    def chain(self, x):

        # Note reverse order compared to AMatrix.
        return x * self


    def cascade(self, x):

        return self.chain(x)


class GMatrix(_TwoPortMatrix):
    """
    +-  -+     +-       -+   +-  -+
    | V2 |  =  | G11  G12|   | I2 |
    | I1 |     | G21  G22|   | V1 |
    +-  -+     +-       -+   +-  -+

           +-         -+
    units  | ohm     1 |
           | 1   1/ohm |
           +-         -+ 

    G = inv(H)
    """

    @property
    def G11(self):
        return self[0, 0]


    @property
    def G12(self):
        return self[0, 1]


    @property
    def G21(self):
        return self[1, 0]


    @property
    def G22(self):
        return self[1, 1]


    @property
    def A(self):
        #return self.H.A
        return AMatrix(1 / self.G21, self.G22 / self.G21, self.G11 / self.G21, self.det() / self.G21)


    @property
    def B(self):
        #return self.H.B
        return AMatrix(-self.det() / self.G12, self.G22 / self.G12, self.G11 / self.G12, -1 / self.G12)


    @property
    def G(self):
        # Perhaps we should make a copy?
        return self


    @property
    def H(self):
        return HMatrix(self.inv())


    @property
    def Y(self):
        return self.H.Y


    @property
    def Z(self):
        return self.H.Z


class HMatrix(_TwoPortMatrix):
    """
    +-  -+     +-       -+   +-  -+
    | V1 |  =  | H11  H12|   | I1 |
    | I2 |     | H21  H22|   | V2 |
    +-  -+     +-       -+   +-  -+

           +-         -+
    units  | ohm     1 |
           | 1   1/ohm |
           +-         -+ 

    H = inv(G)
    """

    @property
    def H11(self):
        return self[0, 0]


    @property
    def H12(self):
        return self[0, 1]


    @property
    def H21(self):
        return self[1, 0]


    @property
    def H22(self):
        return self[1, 1]


    @property
    def A(self):
        return AMatrix(-self.det() / self.H21, -self.H11 / self.H21, -self.H22 / self.H21, -1 / self.H21)


    @property
    def B(self):
        return BMatrix(1 / self.H12, -self.H11 / self.H12, -self.H22 / self.H12, self.det() / self.H12)


    @property
    def H(self):
        # Perhaps we should make a copy?
        return self


    @property
    def Y(self):
        return YMatrix(1 / self.H11, -self.H12 / self.H11, self.H21 / self.H11, self.det() / self.H11)


    @property
    def Z(self):
        return ZMatrix(self.det() / self.H22, self.H12 / self.H22, -self.H21 / self.H22, 1 / self.H22)


class YMatrix(_TwoPortMatrix):
    """
    +-  -+     +-       -+   +-  -+
    | I1 |  =  | Y11  Y12|   | V1 |
    | I2 |     | Y21  Y22|   | V2 |
    +-  -+     +-       -+   +-  -+

           +-           -+
    units  | 1/ohm 1/ohm |
           | 1/ohm 1/ohm |
           +-           -+ 

    Y = inv(Z)
    """

    @property
    def Y11(self):
        return self[0, 0]


    @property
    def Y12(self):
        return self[0, 1]


    @property
    def Y21(self):
        return self[1, 0]


    @property
    def Y22(self):
        return self[1, 1]


    @property
    def Ysc(self):
        return YVector(self.Y[0, 0], self.Y[1, 1])


    @property
    def A(self):
        return AMatrix(-self.Y22 / self.Y21, -1 / self.Y21, -self.det() / self.Y21, -self.Y11 / self.Y21)


    @property
    def B(self):
        return BMatrix(-self.Y11 / self.Y12, 1 / self.Y12, self.det() / self.Y12, -self.Y22 / self.Y12)


    @property
    def H(self):
        return HMatrix(1 / self.Y11, -self.Y12 / self.Y11, self.Y21 / self.Y11, self.det() / self.Y11)


    @property
    def Y(self):
        # Perhaps we should make a copy?
        return self


    @property
    def Z(self):
        # Inverse
        det = self.det()
        return ZMatrix(self.Y22 / det, -self.Y12 / det, -self.Y21 / det, self.Y11 / det)


class ZMatrix(_TwoPortMatrix):
    """
    +-  -+     +-       -+   +-  -+
    | V1 |  =  | Z11  Z12|   | I1 |
    | V2 |     | Z21  Z22|   | I2 |
    +-  -+     +-       -+   +-  -+

           +-         -+
    units  | ohm   ohm |
           | ohm   ohm |
           +-         -+ 

    Z = inv(Y)
    """

    @property
    def Z11(self):
        return self[0, 0]


    @property
    def Z12(self):
        return self[0, 1]


    @property
    def Z21(self):
        return self[1, 0]


    @property
    def Z22(self):
        return self[1, 1]


    @property
    def Zoc(self):
        return ZVector(self.Z[0, 0], self.Z[1, 1])


    @property
    def A(self):
        return AMatrix(self.Z11 / self.Z21, self.det() / self.Z21, 1 / self.Z21, self.Z22 / self.Z21)


    @property
    def B(self):
        return BMatrix(self.Z22 / self.Z12, -self.det() / self.Z12, -1 / self.Z12, self.Z11 / self.Z12)


    @property
    def H(self):
        return HMatrix(self.det() / self.Z22, self.Z12 / self.Z22, -self.Z21 / self.Z22, 1 / self.Z22)


    @property
    def Y(self):
        # Inverse
        det = self.det()
        return YMatrix(self.Z22 / det, -self.Z12 / det, -self.Z21 / det, self.Z11 / det)


    @property
    def Z(self):
        # Perhaps we should make a copy?
        return self


    @classmethod
    def Lsection(cls, Z1, Z2):
        return cls.Tsection(Z1 + Z2, Z2, Z2, Z2)


    @classmethod
    def Tsection(cls, Z1, Z2, Z3):
        # Note if Z3 is infinity then all elements of Z are infinite.
        # Thus we cannot model a single series R with a Z matrix.
        # A single shunt R works though.
        return cls(Z1 + Z2, Z2, Z2, Z2 + Z3)
    

    @classmethod
    def Pisection(cls, Z1, Z2, Z3):

        Za, Zb, Zc = DeltaWye(Z1, Z2, Z3)
        return cls.Tsection(Za, Zb, Zc)



class TwoPort(object):
    """
    two-port networks are constrained to have the same current at each
    port (but flowing in opposite directions).  This is called the
    port condition.

    """

    @property
    def isbuffered(self):
        """Return true if two-port is buffered, i.e., any load
        on the output has no affect on the input. """
        #return self.A.A12 == 0 and self.A.A22 == 0
        return self.B.B12 == 0 and self.B.B22 == 0


    @property
    def isbilateral(self):
        """Return true if two-port is bilateral. """
        return self.B.det() == 1


    @property
    def issymmetrical(self):
        """Return true if two-port is symmetrical. """
        return self.B.B11 == self.B.B22


    @property
    def isseries(self):
        """Return true if two-port is a series network. """
        #return (self.A.A11 == 1) and (self.A.A22 == 1) and (self.A.A21 == 0)
        return (self.B.B11 == 1) and (self.B.B22 == 1) and (self.B.B21 == 0)


    @property
    def isshunt(self):
        """Return true if two-port is a shunt network. """
        #return (self.A.A11 == 1) and (self.A.A22 == 1) and (self.A.A12 == 0)
        return (self.B.B11 == 1) and (self.B.B22 == 1) and (self.B.B12 == 0)


    @property
    def A(self):    
        """Return chain matrix"""
        return self._M.A


    @property
    def B(self):    
        """Return inverse chain matrix"""
        return self._M.B


    @property
    def G(self):    
        """Return inverse hybrid matrix"""
        return self._M.G


    @property
    def H(self):    
        """Return hybrid matrix"""
        return self._M.H


    @property
    def Y(self):    
        """Return admittance matrix"""
        return self._M.Y


    @property
    def Z(self):    
        """Return impedance matrix"""
        return self._M.Z


    @property
    def V1h(self):    
        return Vs(-self.V2b / self.B.B11)


    @property
    def I2h(self):    
        return Is(-self.V2b * self.B.B21 / self.B.B11) - self.I2b


    @property
    def I1y(self):    
        return Is(-self.V2b / self.B.B12)


    @property
    def I2y(self):    
        return Is(self.V2b * self.B.B22 / self.B.B12) - self.I2b


    @property
    def V1z(self):    
        return Vs(-self.I2b / self.B.B21)


    @property
    def V2z(self):    
        return self.V2b - Vs(self.I2b * self.B.B11 / self.B.B21)


    @property
    def Yoc(self):    
        """Return admittance vector with ports open circuit"""
        return YVector(Ys(1 / self.Z1oc), Ys(1 / self.Z2oc))


    @property
    def Y1oc(self):    
        """Return input impedance with the output port open circuit"""
        return Zs(1 / self.Z1oc)


    @property
    def Y2oc(self):    
        """Return output impedance with the input port open circuit"""
        return Ys(1 / self.Z2oc)


    @property
    def Ysc(self):    
        """Return admittance vector with ports short circuit"""
        return self.Y.Ysc


    @property
    def Y1sc(self):    
        """Return input admittance with output port short circuit"""
        return Ys(self.Ysc[0])


    @property
    def Y2sc(self):    
        """Return output admittance with output port short circuit"""
        return Ys(self.Ysc[1])


    @property
    def Zoc(self):    
        """Return impedance vector with ports open circuit"""
        return self.Z.Zoc


    @property
    def Z1oc(self):    
        """Return input impedance with the output port open circuit"""
        return Zs(self.Zoc[0])


    @property
    def Z2oc(self):    
        """Return output impedance with the input port open circuit"""
        return Zs(self.Zoc[1])


    @property
    def Zsc(self):    
        """Return impedance vector with ports short circuit"""
        return ZVector(Zs(1 / self.Y1sc), Zs(1 / self.Y2sc))


    @property
    def Z1sc(self):    
        """Return input impedance with the output port short circuit"""
        return Zs(1 / self.Y1sc)


    @property
    def Z2sc(self):    
        """Return output impedance with the input port short circuit"""
        return Zs(1 / self.Y2sc)


    @property
    def Vgain12(self):
        """Return V2 / V1 for I2 = 0 (forward voltage gain) with
        internal sources zero

        Av = G21 = 1 / A11 = -|B| / B22 = Z21 / Z11 =  Y21 / Y22
        """

        return self.Vgain(1, 2)


    @property
    def Vtransfer(self):
        """Return V2 / V1 for I2 = 0 (forward voltage gain) with
        internal sources zero  (see Vgain12)"""

        return self.Vgain12


    @property
    def Igain12(self):
        """Return I2 / I1 for V2 = 0 (forward current gain) with
        internal sources zero

        Ai = H21 = -1 / A22 = -|B| / B11 = -Z21 / Z22 = Y21 / Y11
        """

        return self.Igain(1, 2)


    @property
    def Itransfer(self):
        """Return I2 / I1 for V2 = 0 (forward current gain) with
        internal sources zero  (sett Igain12)"""

        return self.Igain12


    def Vresponse(self, V, inport=1, outport=2):
        """Return voltage response for specified applied voltage and
        specified ports"""

        if issubclass(V.__class__, OnePort):
            V = V.V
        
        p1 = inport - 1
        p2 = outport - 1

        return Vs(self.Voc[p2] + (V - self.Voc[p1]) * self.Z[p2, p1] / self.Z[p1, p1])


    def Iresponse(self, I, inport=1, outport=2):
        """Return current response for specified applied current and
        specified ports"""

        if issubclass(I.__class__, OnePort):
            I = I.I

        p1 = inport - 1
        p2 = outport - 1

        Y = self.Y
        Isc = self.Isc
                
        return Is(Isc[p2] + Y[p2, p1] / Y[p1, p1] * (I - Isc[p1]))


    def Ytrans(self, inport=1, outport=2):
        """Return transadmittance for specified ports with internal
        sources zero"""

        return Ys(self.Y[outport - 1, inport - 1])


    @property
    def Ytrans12(self):
        """Return I2 / V1 for V2 = 0 (forward transadmittance) with
         internal sources zero

         Y21 = -1 / A12 = |B| / B12
         """

        return Ys(self.Y.Y21)


    @property
    def Ytransfer(self):
        """Return I2 / V1 for V2 = 0 (forward transadmittance) with
         internal sources zero.  This is an alias for Ytrans12.

         Y21 = -1 / A12 = |B| / B12
         """

        return self.Ytrans12


    def Ztrans(self, inport=1, outport=2):
        """Return transimpedance for specified ports with internal
        sources zero"""

        return Zs(self.Z[outport - 1, inport - 1])


    def Ztrans12(self):
        """Return V2 / I1 for I2 = 0 (forward transimpedance) with
        internal sources zero

        Z21 = 1 / A21 = -|B| / B21
        """

        return Zs(self.Z.Z21)


    @property
    def Ztransfer(self):
        """Return V2 / I1 for I2 = 0 (forward transimpedance) with
        internal sources zero.  This is an alias for Ztrans12.

        Z21 = 1 / A21 = -|B| / B21
        """

        return self.Ztrans12


    @property
    def V1oc(self):    
        """Return V1 with all ports open-circuited (i.e., I1 = I2 = 0)"""
        return Vs(self.Voc[0])


    @property
    def V2oc(self):    
        """Return V2 with all ports open-circuited (i.e., I1 = I2 = 0)"""
        return Vs(self.Voc[1])


    @property
    def I1sc(self):    
        """Return I1 with all ports short-circuited, i.e, V1 = V2 = 0"""
        return Is(self.Isc[0])


    @property
    def I2sc(self):    
        """Return I2 with all ports short-circuited, i.e, V1 = V2 = 0"""
        return Is(self.Isc[1])


    @property
    def Voc(self):    
        """Return voltage vector with all ports open-circuited (i.e., In = 0)"""
        return VsVector(self.V1z, self.V2z)


    @property
    def Isc(self):    
        """Return current vector with all ports short-circuited (i.e., V1 = V2 = 0)"""
        return IsVector(self.I1y, self.I2y)


    def __str__(self):

        return self.A.__str__() + '\n\n' + _strpair('V2b', self.V2b, 'I2b', self.I2b)


    @property
    def Bmodel(self):

        return TwoPortBModel(self.B, self.V2b, self.I2b)


    @property
    def Hmodel(self):

        return TwoPortHModel(self.H, self.V1h, self.I2h)


    @property
    def Ymodel(self):

        if self.isshunt:
            warn('Converting a shunt two-port to a Y model is dodgy...')
        return TwoPortYModel(self.Y, self.I1y, self.I2y)


    @property
    def Zmodel(self):

        if self.isseries:
            warn('Converting a series two-port to a Z model is dodgy...')
        return TwoPortZModel(self.Z, self.V1z, self.V2z)


    def chain(self, x):
        """Return the model with, x, appended (cascade or chain connection)"""

        if not issubclass(x.__class__, TwoPort):
            raise TypeError('Argument not', TwoPort)

        # The voltage and current sources can be transformed from the
        # input of a network to its output using:
        # 
        # +-   +     +-       -+   +-  -+
        # | V2b |  =  | B11  B12|   |-Vi |
        # | I2b |     | B21  B22|   | Ii |
        # +-  -+     +-       -+   +-  -+            
        #
        # where the positive connection of Vi is connected to the input pin.

        B = x.B
        foo = B * np.matrix((-self.V2b, self.I2b)).T

        B2 = B * self.B
        return TwoPortBModel(B2, V2b=Vs(foo[0, 0]) + x.V2b, I2b=Is(foo[1, 0]) + x.I2b)


    def append(self, x):
        """Return the model with, x, appended"""

        return self.chain(x)


    def prepend(self, x):
        """Return the model with, x, prepended"""

        return x.chain(self)


    def cascade(self, x):
        """Return the model with, x, appended"""

        return self.chain(x)


    def series(self, x, port=None):
        """Return the model with, x, in series.

         In general, this is tricky to ensure that the port condition
         is valid.  The common ground connection of the first two-port
         shorts out the top of the T of the second two-port.
         """

        if issubclass(x.__class__, OnePort):
            raise NotImplementedError('TODO')

        if issubclass(x.__class__, OnePort):
            raise NotImplementedError('TODO')

        if not issubclass(x.__class__, TwoPort):
            raise TypeError('Argument not ', TwoPort)

        warn('Do you mean chain?  The result of a series combination of two two-ports may be dodgy')
        
        V1z = self.V1z + x.V1z
        V2z = self.V2z + x.V2z
        Z = self.Z + x.Z
        return TwoPortZModel(Z, V1z, V2z)


    def terminate(self, T, port=2):
        """Connect one-port in parallel to specified port and return a Thevenin
        (one-port) object"""

        if port == 1:
            return self.source(T)
        if port == 2:
            return self.load(T)
        raise ValueError('Invalid port ' + port)


    def parallel(self, x, port=None):
        """Return the model with, x, in parallel"""

        if issubclass(x.__class__, OnePort):
            raise NotImplementedError('TODO')

        if not issubclass(x.__class__, TwoPort):
            raise TypeError('Argument not ', TwoPort)
        
        I1y = self.I1y + x.I1y
        I2y = self.I2y + x.I2y
        Y = self.Y + x.Y
        return TwoPortYModel(Y, I1y, I2y)


    def hybrid(self, x, port=None):
        """Return the model with, x, in hybrid connection (series
        input, parallel output)"""

        if issubclass(x.__class__, OnePort):
            raise NotImplementedError('TODO')

        if not issubclass(x.__class__, TwoPort):
            raise TypeError('Argument not ', TwoPort)
        
        V1h = self.V1h + x.V1h
        I2h = self.I2h + x.I2h
        H = self.H + x.H
        return TwoPortHModel(H, V1h, I2h)


    def inversehybrid(self, x, port=None):
        """Return the model with, x, in inverse hybrid connection
        (parallel input, series output)"""

        if issubclass(x.__class__, OnePort):
            raise NotImplementedError('TODO')

        if not issubclass(x.__class__, TwoPort):
            raise TypeError('Argument not ', TwoPort)
        
        raise NotImplementedError('TODO GModel')

        V1g = self.V1g + x.V1g
        I2g = self.I2g + x.I2g
        G = self.G + x.G
        return TwoPortGModel(G, V1g, I2g)


    # Other operations: swapping the input terminals negates the A matrix.
    # switching ports.


    def bridge(self, T):
        """Bridge the ports with a one-port element"""

        if not issubclass(T.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        # FIXME
        return self.parallel(Series(T))


    def load(self, T):
        """Apply a one-port load and return a Thevenin (one-port) object"""

        if not issubclass(T.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        foo = self.chain(Shunt(T))
        return Thevenin(Zs(foo.Z1oc), foo.V1oc)


    def source(self, T):
        """Apply a one-port source and return a Thevenin (one-port) object"""

        if not issubclass(T.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        foo = Shunt(T).chain(self)
        return Thevenin(Zs(foo.Z2oc), foo.V2oc)


    def short_circuit(self, port=2):
        """Apply a short-circuit to specified port and return a
        one-port object"""

        p = port - 1
        Y = self.Y[1 - p, 1 - p]
        I = self.Isc[1 - p]

        return Norton(Ys(Y), Is(I))


    def open_circuit(self, port=2):
        """Apply a open-circuit to specified port and return a
        one-port object"""

        p = port - 1
        Z = self.Z[1 - p, 1 - p]
        V = self.Voc[1 - p]

        return Thevenin(Zs(Z), Vs(V))


class TwoPortBModel(TwoPort):
    """
            +-------------------+    +------+
     I1     |                   | I2'| -  + |          I2
    -->-----+                   +-<--+  V2b +----+-----<--
            |    two-port       |    |      |    |
    +       |    network        | +  +------+    |       +
    V1      |    without        | V2'        +---+---+  V2
    -       |    sources        | -          |   |   |   -
            |    represented    |               I2b  |    
            |    by B matrix    |            |   v   |
            |                   |            +---+---+
            |                   |                |    
    --------+                   +----------------+--------
            |                   |
            +-------------------+

    +-   +     +-        -+   +-  -+     +-   -+
    | V2 |  =  | B11  B12 |   | V1 |  +  | V2b |
    |-I2 |     | B21  B22 |   | I1 |     | I2b |
    +-  -+     +-        -+   +-  -+     +-   -+ 

    
    +-    +     +-        -+   +-  -+
    | V2' |  =  | B11  B12 |   | V1 |
    |-I2' |     | B21  B22 |   | I1 |
    +-   -+     +-        -+   +-  -+    

    +-    +     +-  -+    +-   -+
    | V2' |  =  | V2 | -  | V2b |
    | I2' |     | I2'|    | I2b |
    +-  - +     +-  -+    +-   -+    

    +-         +     +-        -+   +-  -+
    | V2 - V2b |  =  | B11  B12 |   | V1 |
    |-I2 - I2b |     | B21  B22 |   | I1 |
    +-        -+     +-        -+   +-  -+    

    +-   +     +-        -+   +-    +
    | V1 |  =  | A11  A12 |   | V2' |
    | I1 |     | A21  A22 |   |-I2' |
    +-  -+     +-        -+   +-   -+    

    +-   +     +-        -+   +-       -+
    | V1 |  =  | A11  A12 |   | V2 - V2b |
    | I1 |     | A21  A22 |   |-I2 + I2b |
    +-  -+     +-        -+   +-        -+    

    """

    # A disadvantage of the Y and Z matrices is that they become
    # singular for some simple networks.  For example, the Z matrix is
    # singular for a shunt element and the Y matrix is singular for a
    # series element.  The A and B matrices do not seem to have this
    # problem, however, they cannot be extended to three or more ports.
    #

    def __init__(self, B, V2b=Vs(0), I2b=Is(0)):

        if issubclass(B.__class__, TwoPortBModel):
            B, V2b, I2b = B._M, B.Vs2b, B.Is2b, 

        if not isinstance(B, BMatrix):
            raise ValueError('B not BMatrix')

        if not isinstance(V2b, Vs):
            raise ValueError('V2b not Vs')

        if not isinstance(I2b, Is):
            raise ValueError('I2b not Is')
        
        self._M = B
        self.Vs2b = V2b
        self.Is2b = I2b


    def __str__(self):

        return self.B.__str__() + '\n\n' + _strpair('V2b', self.V2b, 'I2b', self.I2b)


    @property
    def B(self):    
        """Return chain matrix"""
        return self._M


    @property
    def I2b(self):    
        return self.Is2b


    @property
    def V2b(self):    
        return self.Vs2b


    @property
    def V1h(self):    
        return Vs(-self.V2b / self.B.B11)


    @property
    def I2h(self):    
        return Is(-self.V2b * self.B.B21 / self.B.B11) - self.I2b


    @property
    def I1y(self):    
        return Is(-self.V2b / self.B.B12)


    @property
    def I2y(self):    
        return Is(self.V2b * self.B.B22 / self.B.B12) - self.I2b


    @property
    def V1z(self):    
        return Vs(-self.I2b / self.B.B21)


    @property
    def V2z(self):    
        return self.V2b - Vs(self.I2b * self.B.B11 / self.B.B21)


    def Vgain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        # Av  = G21 = 1 / A11 = -|B| / B22 = Z21 / Z11 =  Y21 / Y22
        # Av' = H12 = 1 / B11 =  |A| / A22 = Z12 / Z22 = -Y12 / Y11

        if inport == outport:
            return Avs(1)
        if inport == 1 and outport == 2:
            return Avs(1 / self.A.A11)
        if inport == 2 and outport == 1:
            return Avs(1 / self.B.B11)
        raise ValueError('bad port values')


    def Igain(self, inport=1, outport=2):
        """Return current gain for specified ports with internal
         sources zero"""

        # Ai  = H21 = -1 / A22 = -|B| / B11 = -Z21 / Z22 = Y21 / Y11
        # Ai' = G12 =  1 / B22 =  |A| / A11 = -Z12 / Z11 = Y12 / Y22

        if inport == outport:
            return Ais(1)
        if inport == 1 and outport == 2:
            return Ais(-1 / self.A.A22)
        if inport == 2 and outport == 1:
            return Ais(-1 / self.B.B22)
        raise ValueError('bad port values')


class TwoPortHModel(TwoPort):
    """
         +------+   +-------------------+    
     I1  | +  - |   |                   | I2'          I2
    -->--+  V1h +---+                   +-<-------+-----<--
         |      |   |    two-port       |         |
    +    +------+   |    network        | +       |       +
    V1              |    without        | V2' ---+---+  V2
    -               |    sources        | -  |   |   |   -
                    |    represented    |    |  I2h  |    
                    |    by H matrix    |    |   v   |
                    |                   |    +---+---+
                    |                   |        |    
    ----------------+                   +--------+--------
                    |                   |
                    +-------------------+


    +-   +     +-        -+   +-  -+     +-   -+
    | V1 |  =  | H11  H12 |   | I1 |  +  | V1h |
    | I2 |     | H21  H22 |   | V2 |     | I2h |
    +-  -+     +-        -+   +-  -+     +-   -+ 
    """

    def __init__(self, H, V1h=Vs(0), I2h=Is(0)):

        if issubclass(H.__class__, TwoPortHModel):
            H, V1h, I2h = H._M, H.Vs1h, H.Is2h, 

        if not isinstance(H, HMatrix):
            raise ValueError('H not HMatrix')

        if not isinstance(V1h, Vs):
            raise ValueError('V1h not Vs')

        if not isinstance(I2h, Is):
            raise ValueError('I2h not Is')
        
        self._M = H
        self.Vs1h = V1h
        self.Is2h = I2h


    def __str__(self):

        return self.H.__str__() + '\n\n' + _strpair('V1h', self.V1h, 'I2h', self.I2h)


    @property
    def H(self):    
        """Return hybrid matrix"""
        return self._M


    @property
    def V2b(self):    
        """Return V2b"""

        return Vs(self.V1h / self.H.H12)


    @property
    def I2b(self):    
        """Return I2b"""

        return Is(self.H.H22 / self.H.H12 * self.V1h) - self.I2h


    @property
    def V1h(self):    
        return self.Vs1h


    @property
    def I2h(self):    
        return self.Is2h


    def Vgain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        # Av  = G21 = 1 / A11 = -|B| / B22 = Z21 / Z11 =  Y21/Y22
        # Av' = H12 = 1 / B11 =  |A| / A22 = Z12 / Z22 = -Y12/Y11

        if inport == outport:
            return Avs(1)
        if inport == 1 and outport == 2:
            return Avs(self.G.G21)
        if inport == 2 and outport == 1:
            return Avs(self.H.H12)
        raise ValueError('bad port values')


    def Igain(self, inport=1, outport=2):
        """Return current gain for specified ports with internal
        sources zero"""

        # Ai  = H21 = -1 / A22 = -|B| / B11 = -Z21 / Z22 = Y21/Y11
        # Ai' = G12 =  1 / B22 =  |A| / A11 = -Z12 / Z11 = Y12/Y22

        if inport == outport:
            return Ais(1)
        if inport == 1 and outport == 2:
            return Ais(self.H.H21)
        if inport == 2 and outport == 1:
            return Ais(self.G.G12)
        raise ValueError('bad port values')

    
class TwoPortYModel(TwoPort):
    """
                     +-------------------+ 
     I1              |                   | I2'           I2
    -->----+---------+                   +-<-------+-----<--
           |         |    two-port       |         |
    +      |       + |    network        | +       +       +
    V1 +---+---+  V1'|    without        | V2' +---+---+  V2
    -  |   |   |   - |    sources        | -   |   |   |   -
       |  I1y  |     |    represented    |     |  I2y  |    
       |   v   |     |    by Y matrix    |     |   v   |
       +---+---+     |                   |     +---+---+
           |         |                   |         |    
    -------+---------+                   +---------+--------
                     |                   |
                     +-------------------+

    +-   +     +-        -+   +-  -+     +-   -+
    | I1 |  =  | Y11  Y12 |   | V1 |  +  | I1y |
    | I2 |     | Y21  Y22 |   | V2 |     | I2y |
    +-  -+     +-        -+   +-  -+     +-   -+ 

    Ymn = Im / Vn for Vm = 0

    """

    def __init__(self, Y, I1y=Is(0), I2y=Is(0)):

        if issubclass(Y.__class__, TwoPortYModel):
            Y, I1y, I2y = Y._M, Y.Is1y, Y.Is2y

        if not isinstance(Y, YMatrix):
            raise ValueError('Y not YMatrix')

        if not isinstance(I1y, Is):
            raise ValueError('I1y not Is')
        if not isinstance(I2y, Is):
            raise ValueError('I2y not Is')
        
        self._M = Y
        self.Is1y = I1y
        self.Is2y = I2y


    def __str__(self):

        return self.Y.__str__() + '\n\n' + _strpair('I1y', self.Is1y, 'I2y', self.Is2y)


    @property
    def Y(self):    
        """Return admittance matrix"""
        return self._M


    @property
    def I2b(self):    
        return Is(-self.I1y * self.Y.Y11 * self.Y.Y22 / self.Y.Y12) - self.I2y


    @property
    def V2b(self):    
        return Vs(self.I1y * self.Y.Y11 / self.Y.Y22)


    @property
    def I1y(self):    
        return self.Is1y


    @property
    def I2y(self):    
        return self.Is2y

        

class TwoPortZModel(TwoPort):
    """

         +------+    +-------------------+    +------+
    I1   | +  - | I1'|                   | I2'| -  + |  I2
    -->--+  V1z +-->-+                   +-<--+  V2z +--<--
         |      |    |    two-port       |    |      | 
    +    +------+  + |    network        | +  +------+    +
    V1            V1'|    without        | V2'           V2
    -              - |    sources        | -              -
                     |    represented    |             
                     |    by Z matrix    |           
                     |                   |           
                     |                   |          
    -----------------+                   +-----------------
                     |                   |
                     +-------------------+

    +-   +     +-        -+   +-  -+     +-   -+
    | V1 |  =  | Z11  Z12 |   | I1 |  +  | V1z |
    | V2 |     | Z21  Z22 |   | I2 |     | V2z |
    +-  -+     +-        -+   +-  -+     +-   -+ 

    """
        
    def __init__(self, Z, V1z=Vs(0), V2z=Vs(0)):

        if issubclass(Z.__class__, TwoPortZModel):
            Z, V1z, V2z = Z._M, Z.Vs1z, Z.Vs2z

        if not isinstance(Z, ZMatrix):
            raise ValueError('Z not ZMatrix')

        if not isinstance(V1z, Vs):
            raise ValueError('V1z not Vs')
        if not isinstance(V2z, Vs):
            raise ValueError('V2z not Vs')
        
        self._M = Z
        self.Vs1z = V1z
        self.Vs2z = V2z


    def __str__(self):

        return self.Z.__str__() + '\n\n' + _strpair('V1z', self.Vs1z, 'V2z', self.Vs2z)


    @property
    def Z(self):    
        """Return impedance matrix"""
        return self._M


    @property
    def I2b(self):    
        return Is(self.V1z / self.Z.Z12)


    @property
    def V2b(self):    
        return self.V2z - Vs(self.V1z * self.Z.Z22 / self.Z.Z12)


    @property
    def I1y(self):    
        
        Zdet = self.Z.det()
        return Is(-self.V1z * self.Z.Z22 / Zdet - self.V2z * self.Z.Z12 / Zdet)


    @property
    def I2y(self):    

        Zdet = self.Z.det()
        return Is(self.V1z * self.Z.Z21 / Zdet - self.V2z * self.Z.Z11 / Zdet)


    @property
    def V1z(self):    
        return self.Vs1z


    @property
    def V2z(self):    
        return self.Vs2z


    def Vgain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        # Av  = G21 = 1 / A11 = -|B| / B22 = Z21 / Z11 =  Y21/Y22
        # Av' = H12 = 1 / B11 =  |A| / A22 = Z12 / Z22 = -Y12/Y11

        if inport == outport:
            return Avs(1)
        if inport == 1 and outport == 2:
            return Avs(self.Z.Z21 / self.Z.Z11)
        if inport == 2 and outport == 1:
            return Avs(self.Z.Z12 / self.Z.Z22)
        raise ValueError('bad port values')


    def Igain(self, inport=1, outport=2):
        """Return current gain for specified ports with internal
        sources zero"""

        # Ai  = H21 = -1 / A22 = -|B| / B11 = -Z21 / Z22 = Y21/Y11
        # Ai' = G12 =  1 / B22 =  |A| / A11 = -Z12 / Z11 = Y12/Y22

        if inport == outport:
            return Ais(1)
        if inport == 1 and outport == 2:
            return Ais(-self.Z.Z21 / self.Z.Z22)
        if inport == 2 and outport == 1:
            return Ais(-self.Z.Z12 / self.Z.Z11)
        raise ValueError('bad port values')


class Series(TwoPortBModel):

    def __init__(self, OP1):
        """
           +---------+   
         --+   OP1   +---
           +---------+   

         ----------------
         """

        super (Series, self).__init__(BMatrix.Zseries(OP1.Z), V2b=OP1.V)


class Shunt(TwoPortBModel):

    def __init__(self, OP1):
        """
                 
         -----+----
              |    
            +-+-+  
            |   |  
            |OP1|  
            |   |  
            +-+-+  
              |    
         -----+----
         """

        super (Shunt, self).__init__(BMatrix.Yshunt(OP1.Y), V2b=Vs(0), I2b=Is(OP1.I))


class IdealTransformer(TwoPortBModel):
    """Ideal transformer voltage gain alpha, current gain 1 / alpha"""


    def __init__(self, alpha=1):


        super (IdealTransformer, self).__init__(BMatrix.transformer(alpha))


class IdealGyrator(TwoPortBModel):
    """Ideal gyrator with gyration resistance R.
    
    A gyrator converts a voltage to current and a current to voltage.
    Cascaded gyrators act like a transformer"""


    def __init__(self, R=1):


        super (IdealGyrator, self).__init__(BMatrix.gyrator(R))


class VoltageFollower(TwoPortBModel):
    """Voltage follower"""

    def __init__(self):

        super (VoltageFollower, self).__init__(BMatrix.voltage_amplifier(1))


class VoltageAmplifier(TwoPortBModel):
    """Voltage amplifier"""

    def __init__(self, Av=1, Af=0, Yin=0, Zout=0):

        super (VoltageAmplifier, self).__init__(BMatrix.voltage_amplifier(Av, Af, Yin, Zout))


class IdealVoltageAmplifier(TwoPortBModel):
    """Ideal voltage amplifier"""

    def __init__(self, Av=1):

        super (IdealVoltageAmplifier, self).__init__(BMatrix.voltage_amplifier(Av))


class IdealDelay(TwoPortBModel):
    """Ideal buffered delay"""

    def __init__(self, delay=0):

        delay = _Expr(delay)
        super (IdealDelay, self).__init__(BMatrix.voltage_amplifier(sym.exp(-delay * delay.s)))


class IdealVoltageDifferentiator(TwoPortBModel):
    """Voltage differentiator"""

    def __init__(self, Av=1):

        super (IdealVoltageDifferentiator, self).__init__(BMatrix.voltage_differentiator(Av))


class IdealVoltageIntegrator(TwoPortBModel):
    """Ideal voltage integrator"""

    def __init__(self, Av=1):

        super (IdealVoltageIntegrator, self).__init__(BMatrix.voltage_integrator(Av))


class CurrentFollower(TwoPortBModel):
    """Current follower"""

    def __init__(self):

        super (CurrentFollower, self).__init__(BMatrix.current_amplifier(1))


class IdealCurrentAmplifier(TwoPortBModel):
    """Ideal Current amplifier"""

    def __init__(self, Ai=1):

        super (IdealCurrentAmplifier, self).__init__(BMatrix.current_amplifier(Ai))


class IdealCurrentDifferentiator(TwoPortBModel):
    """Ideal Current differentiator"""

    def __init__(self, Av=1):

        super (IdealCurrentDifferentiator, self).__init__(BMatrix.current_differentiator(Av))


class IdealCurrentIntegrator(TwoPortBModel):
    """Ideal Current integrator"""

    def __init__(self, Av=1):

        super (IdealCurrentIntegrator, self).__init__(BMatrix.current_integrator(Av))


class OpampInverter(TwoPortBModel):
    """Opamp inverter"""

    def __init__(self, R1, R2):

        # FIXME for initial voltages.
        super (OpampInverter, self).__init__(AMatrix(-R1.Z / R2.Z, 0, -1 / R2.Z, 0).B)


class OpampIntegrator(TwoPortBModel):
    """Inverting opamp integrator"""

    def __init__(self, R1, C1):

        # FIXME for initial voltages.
        super (OpampIntegrator, self).__init__(AMatrix(-R1.Z / C1.Z, 0, -1 / C1.Z, 0).B)


class OpampDifferentiator(TwoPortBModel):
    """Inverting opamp differentiator"""

    def __init__(self, R1, C1):

        # FIXME for initial voltages.
        super (OpampDifferentiator, self).__init__(AMatrix(-R1.Z * C1.Z, 0, -R1.Z, 0).B)


class TSection(TwoPortBModel):

    def __init__(self, OP1, OP2, OP3):
        """
           +---------+       +---------+       
         --+   OP1   +---+---+   OP3   +---
           +---------+   |   +---------+   
                       +-+-+             
                       |   |             
                       |OP2|             
                       |   |             
                       +-+-+             
                         |               
         ----------------+-----------------

         The Z matrix for a resistive T section is
         [ R1 + R2, R2     ]
         [      R2, R2 + R3]

         """

        super (TSection, self).__init__(Series(OP1).chain(Shunt(OP2)).chain(Series(OP3)))



class TwinTSection(TwoPortBModel):

    def __init__(self, OP1a, OP2a, OP3a, OP1b, OP2b, OP3b):
        """
              +---------+       +---------+       
           +--+   OP1a  +---+---+   OP3a  +--+
           |  +---------+   |   +---------+  |
           |              +-+-+              |
           |              |   |              |
           |              |OP2a              |
           |              |   |              |
           |              +-+-+              |
           |                |                |
           |                v                |
           |  +---------+       +---------+  |    
         --+--+   OP1b  +---+---+   OP3b  +--+--
              +---------+   |   +---------+   
                          +-+-+             
                          |   |             
                          |OP2b             
                          |   |             
                          +-+-+             
                            |               
         -------------------+--------------------

         """

        super (TwinTSection, self).__init__(TSection(OP1a, OP2a, OP3a).parallel(TSection(OP1b, OP2b, OP3b)))



class BridgedTSection(TwoPortBModel):

    def __init__(self, OP1, OP2, OP3, OP4):
        """
                       +---------+       
           +-----------+   OP4   +-----------+
           |           +---------+           |
           |                                 |
           |  +---------+       +---------+  |    
         --+--+   OP1b  +---+---+   OP3b  +--+--
              +---------+   |   +---------+   
                          +-+-+             
                          |   |             
                          |OP2b             
                          |   |             
                          +-+-+             
                            |               
         -------------------+--------------------

         """

        super (TwinTSection, self).__init__(TSection(OP1, OP2, OP3).parallel(Series(OP4)))


class PiSection(TwoPortBModel):

    def __init__(self, OP1, OP2, OP3):
        """
                  +---------+       
        -----+----+   OP2    +---+-----
             |    +---------+   |  
           +-+-+              +-+-+
           |   |              |   |
           |OP1|              |OP3|
           |   |              |   |
           +-+-+              +-+-+                       
             |                  |
        -----+------------------+-----

        """

        super (PiSection, self).__init__(Shunt(OP1).chain(Series(OP2)).chain(Shunt(OP3)))



class LSection(TwoPortBModel):

    def __init__(self, OP1, OP2):
        """
           +---------+       
         --+   OP1   +---+----
           +---------+   |   
                       +-+-+ 
                       |   | 
                       |OP2 | 
                       |   | 
                       +-+-+ 
                         |   
         ----------------+----
         """

        super (LSection, self).__init__(Series(OP1).chain(Shunt(OP2)))



class Ladder(TwoPortBModel):

    def __init__(self, OP1, *args):
        """
        
        Create (unbalanced) ladder network with alternating Series and Shunt
        networks chained

           +---------+       +---------+       
         --+   OP1   +---+---+ args[1] +---
           +---------+   |   +---------+   
                       +-+-+             
                       |   |             
                       |   | args[0]             
                       |   |             
                       +-+-+             
                         |               
         ----------------+-----------------
         """

        TP = Series(OP1)

        for m, arg in enumerate(args):
            
            if m & 1:
                TP = TP.chain(Series(arg))
            else:
                TP = TP.chain(Shunt(arg))            

        super (Ladder, self).__init__(TP)


class GeneralTxLine(TwoPortBModel):

    def __init__(self, Z0, gamma, l):
        """
        Z0 is the (real) characteristic impedance (ohms)
        gamma is the propagation constant (1/m)
        l is the transmission line length (m)
        """

        if not _symbolic:
            raise ValueError('Can only implement transmission line symbolically')

        H = sym.exp(gamma * l)

        B11 = 0.5 * (H + 1 / H)
        B12 = 0.5 * (1 / H - H) * Z0
        B21 = 0.5 * (1 / H - H) / Z0
        B22 = 0.5 * (H + 1 / H)

        super (GeneralTxLine, self).__init__(BMatrix(B11, B12, B21, B22))


class LosslessTxLine(GeneralTxLine):

    def __init__(self, Z0, c=1.5e8, l=1):
        """
        Z0 is the (real) characteristic impedance (ohms)
        c is the propagation speed (m/s)
        l is the transmission line length (m)
        """

        if not _symbolic:
            raise ValueError('Can only implement transmission line symbolically')

        s = sym.Symbol('s')
        gamma = s / c

        super (LosslessTxLine, self).__init__(Z0, gamma, l)


class TxLine(GeneralTxLine):

    def __init__(self, R, L, G, C, l=1):
        """
        R series resistance/metre
        L series inductance/metre
        G shunt conductance/metre
        C shunt capacitance/metre
        l is the transmission line length
        """

        if not _symbolic:
            raise ValueError('Can only implement transmission line symbolically')

        s = sym.Symbol('s')

        Z = R + s * L
        Y = G + s * C
        gamma = sym.sqrt(Z * Y)
        Z0 = sym.sqrt(Z / Y)

        super (TxLine, self).__init__(Z0, gamma, l)


class _ThreePortMatrix(sym.Matrix):

    def __new__ (cls, *args):

        if len(args) == 9:
            return super (_ThreePortMatrix, cls).__new__(cls, ((args[0], args[1], args[2]), (args[3], args[4], args[5]), (args[6], args[7], args[8])))

        return super (_ThreePortMatrix, cls).__new__(cls, *args)


    # @property
    # def det(self):
    #     """Return determinant"""


    # def inv(self):
    #     """Return inverse"""


    @property
    def Z(self):
        return ZMatrix3(self.Y.inv())


    @property
    def Y(self):
        return YMatrix3(self.Z.inv())


class ZMatrix3(_ThreePortMatrix):
    """
    +-  -+     +-             -+   +-  -+
    | V1 |     | Z11  Z12  Z13 |   | I1 |
    | V2 |  =  | Z21  Z22  Z23 |   | I2 |
    | V3 |     | Z31  Z32  Z33 |   | I3 |
    +-  -+     +-             -+   +-  -+

    Z = inv(Y)
    """

    @property
    def Z(self):
        # Perhaps we should make a copy?
        return self


class YMatrix3(_ThreePortMatrix):
    """
    +-  -+     +-             -+   +-  -+
    | I1 |  =  | Y11  Y12  Y13 |   | V1 |
    | I2 |     | Y21  Y22  Y23 |   | V2 |
    | I3 |     | Y31  Y32  Y33 |   | V3 |
    +-  -+     +-             -+   +-  -+

    Y = inv(Z)
    """


    @property
    def Y(self):
        # Perhaps we should make a copy?
        return self


class ThreePort(object):
    """ 

    +-  -+     +-             -+   +-  -+     +-   -+
    | V1 |     | Z11  Z12  Z13 |   | I1 |     | V1z |
    | V2 |  =  | Z21  Z22  Z23 |   | I2 |  +  | V2z |
    | V3 |     | Z31  Z32  Z33 |   | I3 |     | V3z |
    +-  -+     +-             -+   +-  -+     +-   -+

    The A, B, G, and H models are invalid for multiports.

    Unfortunately, the Z model can blow up for simple networks.

    """

    def __init__(self, Z, Vz=VsVector((0, 0, 0))):

        if not isinstance(Z, ZMatrix3):
            raise ValueError('Z not ZMatrix3')

        if not isinstance(Vz, VsVector):
            raise ValueError('Vz not VsVector')

        self._M = Z
        self.Vsz = Vz


    @property
    def Voc(self):    
        """Return voltage vector with all ports open-circuited (i.e., In = 0)"""
        return self.Vsz


    @property
    def Isc(self):    
        """Return current vector with all ports short-circuited (i.e., Vn = 0)"""
        Y = self.Y
        Voc = self.Voc

        Isc = IsVector([Voc[m] * Y[m, m] for m in range(len(Voc))])
        return Isc


    @property
    def Y(self):    
        """Return admittance matrix"""
        return YMatrix3(self._M.Y)


    @property
    def Z(self):    
        """Return impedance matrix"""
        return self._M


    @property
    def Yoc(self):    
        """Return admittance vector with ports open circuit"""
        Z = self.Z
        return YVector([1 / Z[m, m] for m in range(Z.shape[0])])


    @property
    def Ysc(self):    
        """Return admittance vector with ports short circuit"""
        Y = self.Y
        return YVector([Y[m, m] for m in range(Y.shape[0])])


    @property
    def Zoc(self):    
        """Return impedance vector with ports open circuit"""
        Z = self.Z
        return ZVector([Z[m, m] for m in range(Z.shape[0])])


    @property
    def Zsc(self):    
        """Return impedance vector with ports short circuit"""
        Y = self.Y
        return ZVector([1 / Y[m, m] for m in range(Y.shape[0])])


    def portcheck(self, port):

        if port not in (1, 2, 3):
            raise ValueError('Invalid port ' + port)


    def Vgain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        self.portcheck(inport)
        self.portcheck(outport)

        p1 = inport - 1
        p2 = outport - 1

        return Avs(self.Z[p2, p1] / self.Z[p1, p1])


    def Igain(self, inport=1, outport=2):
        """Return voltage gain for specified ports with internal
        sources zero"""

        self.portcheck(inport)
        self.portcheck(outport)

        p1 = inport - 1
        p2 = outport - 1

        Y = self.Y

        return Ais(self.Y[p2, p1] / self.Y[p1, p1])


    def Vresponse(self, V, inport=1, outport=2):
        """Return voltage response for specified applied voltage and
        specified ports"""

        self.portcheck(inport)
        self.portcheck(outport)

        p1 = inport - 1
        p2 = outport - 1

        return Vs(self.Voc[p2] + (V - self.Voc[p1]) * self.Z[p2, p1] / self.Z[p1, p1])


    def Iresponse(self, I, inport=1, outport=2):
        """Return current response for specified current voltage and
        specified ports"""

        self.portcheck(inport)
        self.portcheck(outport)

        p1 = inport - 1
        p2 = outport - 1

        Y = self.Y
        Isc = self.Isc
                
        return Is(Isc[p2] + (I - Isc[p1]) * Y[p2, p1] / Y[p1, p1])


    def attach_parallel(self, T, port=2):
        """Attach one-port in parallel to specified port"""

        if not issubclass(T.__class__, OnePort):
            raise TypeError('Argument not ', OnePort)

        self.portcheck(port)

        p = port - 1

        Y = self.Y
        Y[p, p] += T.Y
        Isc = self.Isc
        Isc[p] += T.Isc
        Z = Y.Z
        Voc = VsVector([Vs(Isc[m] * Z[m, m]) for m in range(len(Isc))])
        return ThreePort(Z, Voc)


    def bridge(self, T, inport=1, outport=2):
        """Bridge the specified ports with a one-port element"""

        self.portcheck(inport)
        self.portcheck(outport)

        # Create two-port series element.
        s = Series(T)
        
        # The impedance matrix for a series element is infinite.

        Y3 = YMatrix3(((0, 0, 0), (0, 0, 0), (0, 0, 0)))
        Y2 = s.Y
        p1 = inport - 1
        p2 = outport - 1

        Y3[p1, p1] = Y2[0, 0]
        Y3[p2, p2] = Y2[1, 1]
        Y3[p1, p2] = Y2[0, 1]
        Y3[p2, p1] = Y2[1, 0]

        Y = self.Y + Y3
        Isc = self.Isc
        Isc[p1] -= T.Isc
        Isc[p2] += T.Isc
        Z = Y.Z
        Voc = VsVector([Vs(Isc[m] * Z[m, m]) for m in range(len(Isc))])
        return ThreePort(Y.Z, Voc)

    
    def parallel(self, x, port=None):
        """Return the model with, x, in parallel"""

        if issubclass(x.__class__, OnePort):
            return self.attach_parallel(x, port)

        if issubclass(x.__class__, TwoPort):
            # We could special case a series or shunt network here.
            raise NotImplementedError('TODO')

        if not issubclass(x.__class__, ThreePort):
            raise TypeError('Argument not ', ThreePort)
        
        Y = self.Y + x.Y
        Isc = self.Isc + x.Isc
        Z = Y.Z
        Voc = VsVector([Vs(Isc[m] * Z[m, m]) for m in range(len(Isc))])
        return ThreePort(Z, Voc)


    def series(self, x, port=None):
        """Return the model with, x, in series"""

        if issubclass(x.__class__, OnePort):
            raise NotImplementedError('TODO')

        if issubclass(x.__class__, TwoPort):
            raise NotImplementedError('TODO')

        if not issubclass(x.__class__, ThreePort):
            raise TypeError('Argument not ', ThreePort)
        
        warn('Will this ever work?')

        Z = self.Z + x.Z
        Voc = self.Voc + x.Voc
    
        return ThreePort(Z, Voc)


    def terminate(self, T, port=2):
        """Connect one-port in parallel to specified port and return a
        two-port object"""

        return self.attach_parallel(T, port).open_circuit(port)


    def short_circuit(self, port=2):
        """Apply a short-circuit to specified port and return a
        two-port object"""

        # Remove the unwanted port from the Ymatrix.
        Y = self.Y.copy()
        Y.row_del(port - 1)
        Y.col_del(port - 1)
        Y = YMatrix(Y)

        # CHECKME, perhaps use Isc?
        Voc = self.Voc.copy()
        Voc.row_del(port - 1)

        return TwoPortZModel(Y.Z, Vs(Voc[0]), Vs(Voc[1]))
        

    def open_circuit(self, port=2):
        """Apply a open-circuit to specified port and return a
        two-port object"""

        # Remove the unwanted port from the Zmatrix.        
        Z = self.Z.copy()
        Z.row_del(port - 1)
        Z.col_del(port - 1)
        Z = ZMatrix(Z)

        Voc = self.Voc.copy()
        Voc.row_del(port - 1)

        return TwoPortZModel(Z, Vs(Voc[0]), Vs(Voc[1]))



class Opamp(ThreePort):

    def __init__(self, Rd=1e9, Ro=1e-6, A=100000, Rp=1e9, Rm=1e9):
        """
            |\
            |  \
        1 --+ +  \
            |      \
            |        \
            |          +--- 3
            |        /
            |      /
        2 --+ -  /
            |  /
            |/

        Each port voltage, Vn, is referenced to a common ground.

        """

        # If Ro=0, then Z matrix singular.

        Ra = Rp * (Rd + Rm) / (Rp + Rd + Rm)
        Rb = Rm * (Rd + Rp) / (Rp + Rd + Rm)

        Z = ZMatrix3(((Rp + Rd, Rd, 0),
                     (Rd, Rm + Rd, 0),
                     (A * Ra, -A * Rb, Ro)))
        super (Opamp, self).__init__(Z)
    


def test():

    mC = C(40e-12)
    mL = L(0.06e-4)
    mR = R(2)
    mV = V(5)
    mI = I(20)

    print(mC)
    print(mL)
    print(mR)
    print(mI)
    print(mV)

    mZ1 = mR.Z + mV.V
    mZ2 = mR + mV
    mZ3 = mR | mI

    print(mZ1)
    print(mZ2)
    print(mZ3)

    a = AMatrix.Zseries(Zs(10))
    b = AMatrix.Zshunt(Zs(20))
    c = a.chain(b)

    print(c)

