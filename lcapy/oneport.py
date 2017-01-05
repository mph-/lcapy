"""
This module supports simple linear one-port networks based on the
following ideal components:

V independent voltage source
I independent current source
R resistor
C capacitor
L inductor

These components are converted to s-domain models and so capacitor and
inductor components can be specified with initial voltage and
currents, respectively, to model transient responses.

One-ports can either be connected in series (+) or parallel (|) to
create a new one-port.

Copyright 2014--2017 Michael Hayes, UCECE
"""

from __future__ import division
import sympy as sym
from lcapy.core import t, s, Vs, Is, Zs, Ys, cExpr, sExpr, tExpr, Expr
from lcapy.core import cos, exp, symbol, j, Vphasor, Iphasor, It, Vconst, Iconst, Vn, In
from lcapy.core import Vsuper, Isuper, pretty
from lcapy.sympify import symbols_find
from lcapy.network import Network


__all__ = ('V', 'I', 'v', 'i', 'R', 'L', 'C', 'G', 'Y', 'Z',
           'Vdc', 'Vstep', 'Idc', 'Istep', 'Vac', 'sV', 'sI',
           'Iac', 'Vnoise', 'Inoise', 
           'Par', 'Ser', 'Xtal', 'FerriteBead')

def _check_oneport_args(args):

    for arg1 in args:
        if not isinstance(arg1, OnePort):
            raise ValueError('%s not a OnePort' % arg1)


class OnePort(Network):
    """One-port network

    There are four major types of OnePort:
       VoltageSource
       CurrentSource
       Impedance
       Admittance
       ParSer for combinations of OnePort

    Attributes: Y, Z, Voc, Isc, y, z, voc, isc
      Y = Y(s)  admittance
      Z = Z(s)  impedance
      Voc       open-circuit voltage in appropriate transform domain
      Isc       short-circuit current in appropriate transform domain
      y = y(t)  impulse response of admittance
      z = z(t)  impulse response of impedance
      voc = voc(t) open-circuit voltage time response
      isc = isc(t) short-circuit current time response
    """

    # Dimensions and separations of component with horizontal orientation.
    height = 0.3
    hsep = 0.5
    width = 1
    wsep = 0.5

    netname = ''
    netkeyword = ''

    _Z = None
    _Y = None
    _Voc = None
    _Isc = None

    @property
    def Z(self):
        if self._Z is not None:
            return self._Z
        if self._Y is not None:
            return Zs(1 / self._Y)
        if self._Voc is not None:        
            return Zs(0)
        if self._Isc is not None:        
            return Zs(1 / Ys(0))
        raise ValueError('_Isc, _Voc, _Y, or _Z undefined for %s' % self)

    @property
    def Y(self):
        if self._Y is not None:
            return self._Y
        return Ys(1 / self.Z)

    @property
    def Voc(self):
        if self._Voc is not None:
            return self._Voc
        if self._Isc is not None:
            return self._Isc * self.Z
        if self._Z is not None:        
            return Vsuper(0)
        if self._Y is not None:        
            return Isuper(0)
        raise ValueError('_Isc, _Voc, _Y, or _Z undefined for %s' % self)        

    @property
    def Isc(self):
        if self._Isc is not None:
            return self._Isc
        return self.Voc / self.Z

    def __add__(self, OP):
        """Series combination"""

        return Ser(self, OP)

    def __or__(self, OP):
        """Parallel combination"""

        return Par(self, OP)

    def series(self, OP):
        """Series combination"""

        return Ser(self, OP)

    def parallel(self, OP):
        """Parallel combination"""

        return Par(self, OP)

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

    def expand(self):

        return self

    @property
    def voc(self):
        return self.Voc.time()

    @property
    def isc(self):
        return self.Isc.time()

    @property
    def z(self):
        return self.Z.time()

    @property
    def y(self):
        return self.Y.time()

    def thevenin(self):
        """Simplify to a Thevenin network"""

        Voc = self.Voc
        Z = self.Z
        
        if Voc.is_superposition and not Z.is_real:
            print('Warning, detected superposition with reactive impedance,'
                  ' using s-domain.')
            Voc = Voc.laplace()
        elif Voc.is_ac:
            Z = Z(j * Voc.ac_keys()[0])
        elif Voc.is_dc:
            Z = Z(0)

        V1 = Voc.cpt()
        Z1 = Z.cpt()

        if Voc == 0:
            return Z1
        if Z == 0:
            return V1

        return Ser(Z1, V1)

    def norton(self):
        """Simplify to a Norton network"""

        Isc = self.Isc
        Y = self.Y
        
        if Isc.is_superposition and not Y.is_real:
            print('Warning, detected superposition with reactive impedance,'
                  ' using s-domain.')
            Isc = Isc.laplace()
        elif Isc.is_ac:
            Y = Y(j * Isc.ac_keys()[0])
        elif Isc.is_dc:
            Y = Y(0)

        I1 = Isc.cpt()
        Y1 = Y.cpt()

        if Isc == 0:
            return Y1
        if Y == 0:
            return I1

        return Par(Y1, I1)

    def s_model(self):
        """Convert to s-domain"""

        if self._Voc is not None:
            if self._Voc == 0:
                return Z(self.Z)
            Voc = self._Voc.laplace()
            if self.Z == 0:
                return V(Voc)
            return Ser(V(Voc), Z(self.Z))
        elif self._Isc is not None:
            if self._Isc == 0:
                return Y(self.Y)
            Isc = self._Isc.laplace()
            if self.Y == 0:
                return I(Isc)
            return Par(I(Isc), Y(self.Y))
        elif self._Z is not None:
            return Z(self._Z)
        elif self._Y is not None:
            return Y(self._Y)        
        raise ValueError('Internal error')

        
class ParSer(OnePort):
    """Parallel/serial class"""

    def __str__(self):

        str = ''

        for m, arg in enumerate(self.args):
            argstr = arg.__str__()

            if isinstance(arg, ParSer) and arg.__class__ != self.__class__:
                argstr = '(' + argstr + ')'

            str += argstr

            if m != len(self.args) - 1:
                str += ' %s ' % self._operator

        return str

    def _repr_pretty_(self, p, cycle):

        p.text(self.pretty())

    def _repr_latex_(self):

        return '$%s$' % self.latex()

    def pretty(self):

        str = ''

        for m, arg in enumerate(self.args):
            argstr = arg.pretty()

            if isinstance(arg, ParSer) and arg.__class__ != self.__class__:
                argstr = '(' + argstr + ')'

            str += argstr

            if m != len(self.args) - 1:
                str += ' %s ' % self._operator

        return str

    def pprint(self):

        print(self.pretty())

    def latex(self):

        str = ''

        for m, arg in enumerate(self.args):
            argstr = arg.latex()

            if isinstance(arg, ParSer) and arg.__class__ != self.__class__:
                argstr = '(' + argstr + ')'

            str += argstr

            if m != len(self.args) - 1:
                str += ' %s ' % self._operator

        return str

    def _combine(self, arg1, arg2):

        if arg1.__class__ != arg2.__class__:
            if self.__class__ == Ser:
                if isinstance(arg1, V) and arg1.Voc == 0:
                    return arg2
                if isinstance(arg2, V) and arg2.Voc == 0:
                    return arg1
                if isinstance(arg1, Z) and arg1.Z == 0:
                    return arg2
                if isinstance(arg2, Z) and arg2.Z == 0:
                    return arg1
            if self.__class__ == Par:
                if isinstance(arg1, I) and arg1.Isc == 0:
                    return arg2
                if isinstance(arg2, I) and arg2.Isc == 0:
                    return arg1
                if isinstance(arg1, Y) and arg1.Y == 0:
                    return arg2
                if isinstance(arg2, Y) and arg2.Y == 0:
                    return arg1

            return None

        if self.__class__ == Ser:
            if isinstance(arg1, I):
                return None
            if isinstance(arg1, Vdc):
                return Vdc(arg1.v0 + arg2.v0)
            # Could simplify Vac here if same frequency
            if isinstance(arg1, V):
                return V(arg1 + arg2)
            if isinstance(arg1, R):
                return R(arg1.R + arg2.R)
            if isinstance(arg1, L):
                # The currents should be the same!
                if arg1.i0 != arg2.i0:
                    raise ValueError('Series inductors with different'
                          ' initial currents!')
                return L(arg1.L + arg2.L, arg1.i0)
            if isinstance(arg1, G):
                return G(arg1.G * arg2.G / (arg1.G + arg2.G))
            if isinstance(arg1, C):
                return C(
                    arg1.C * arg2.C / (arg1.C + arg2.C), arg1.v0 + arg2.v0)
            return None

        elif self.__class__ == Par:
            if isinstance(arg1, V):
                return None
            if isinstance(arg1, Idc):
                return Idc(arg1.i0 + arg2.i0)
            # Could simplify Iac here if same frequency
            if isinstance(arg1, I):
                return I(arg1 + arg2)
            if isinstance(arg1, G):
                return G(arg1.G + arg2.G)
            if isinstance(arg1, C):
                # The voltages should be the same!
                if arg1.v0 != arg2.v0:
                    raise ValueError('Parallel capacitors with different'
                          ' initial voltages!')
                return C(arg1.C + arg2.C, arg1.v0)
            if isinstance(arg1, R):
                return R(arg1.R * arg2.R / (arg1.R + arg2.R))
            if isinstance(arg1, L):
                return L(
                    arg1.L * arg2.L / (arg1.L + arg2.L), arg1.i0 + arg2.i0)
            return None

        else:
            raise TypeError('Undefined class')

    def simplify(self, deep=True):
        """Perform simple simplifications, such as parallel resistors,
        series inductors, etc., rather than collapsing to a Thevenin
        or Norton network.

        This does not expand compound components such as crystal
        or ferrite bead models.  Use expand() first.
        """

        # Simplify args (recursively) and combine operators if have
        # Par(Par(A, B), C) etc.
        new = False
        newargs = []
        for m, arg in enumerate(self.args):
            if isinstance(arg, ParSer):
                arg = arg.simplify(deep)
                new = True
                if arg.__class__ == self.__class__:
                    newargs.extend(arg.args)
                else:
                    newargs.append(arg)
            else:
                newargs.append(arg)

        if new:
            self = self.__class__(*newargs)

        # Scan arg list looking for compatible combinations.
        # Could special case the common case of two args.
        new = False
        args = list(self.args)
        for n in range(len(args)):

            arg1 = args[n]
            if arg1 is None:
                continue
            if isinstance(arg1, ParSer):
                continue

            for m in range(n + 1, len(args)):

                arg2 = args[m]
                if arg2 is None:
                    continue
                if isinstance(arg2, ParSer):
                    continue

                # TODO, think how to simplify things such as
                # Par(Ser(V1, R1), Ser(R2, V2)).
                # Could do Thevenin/Norton transformations.

                newarg = self._combine(arg1, arg2)
                if newarg is not None:
                    # print('Combining', arg1, arg2, 'to', newarg)
                    args[m] = None
                    arg1 = newarg
                    new = True

            args[n] = arg1

        if new:
            args = [arg for arg in args if arg is not None]
            if len(args) == 1:
                return args[0]
            self = self.__class__(*args)

        return self

    def expand(self):
        """Expand compound components such as crystals or ferrite bead
        models into R, L, G, C, V, I"""

        newargs = []
        for m, arg in enumerate(self.args):
            newarg = arg.expand()
            newargs.append(newarg)

        return self.__class__(*newargs)

    def s_model(self):
        """Convert to s-domain"""
        args = [arg.s_model() for arg in self.args]
        return (self.__class__(*args))

    @property
    def Isc(self):
        return self.cct.Isc(1, 0)

    @property
    def Voc(self):
        return self.cct.Voc(1, 0)

    @property
    def Y(self):
        # Could extract directly if have Y || I or Z + V
        return self.cct.admittance(1, 0)

    @property
    def Z(self):
        # Could extract directly if have Y || I or Z + V
        return self.cct.impedance(1, 0)

class Par(ParSer):
    """Parallel class"""

    _operator = '|'

    def __init__(self, *args):

        _check_oneport_args(args)
        super(Par, self).__init__()
        self.args = args

        for n, arg1 in enumerate(self.args):
            for arg2 in self.args[n + 1:]:
                if isinstance(arg1, V) and isinstance(arg2, V):
                    raise ValueError('Voltage sources connected in parallel'
                                     ' %s and %s' % (arg1, arg2))
                if isinstance(arg1, V) or isinstance(arg2, V):
                    # Should simplify by removing redundant component to
                    # save later grief with Norton or Thevenin transformation.
                    print('Warn: parallel connection with voltage source:'
                          ' %s and %s' % (arg1, arg2))

    @property
    def width(self):

        total = 0
        for arg in self.args:
            val = arg.width
            if val > total:
                total = val
        return total + 2 * self.wsep

    @property
    def height(self):

        total = 0
        for arg in self.args:
            total += arg.height
        return total + (len(self.args) - 1) * self.hsep

    def net_make(self, net, n1=None, n2=None):

        s = []
        if n1 is None:
            n1 = net.node
        n3, n4 =  net.node, net.node

        H = [(arg.height + self.hsep) * 0.5 for arg in self.args]
        
        N = len(H)
        num_branches = N // 2

        # Draw component in centre if have odd number in parallel.
        if (N & 1):
            s.append(self.args[N // 2].net_make(net, n3, n4))

        na, nb = n3, n4

        s.append('W %s %s; right=%s' % (n1, n3, self.wsep))

        # Draw components above centre
        for n in range(num_branches):

            if not (N & 1) and n == 0:
                sep = H[N // 2 - 1]
            else:
                sep = H[N // 2 - n] + H[N // 2 - 1 - n]

            nc, nd =  net.node, net.node
            s.append('W %s %s; up=%s' % (na, nc, sep))
            s.append('W %s %s; up=%s' % (nb, nd, sep))
            s.append(self.args[N // 2 - 1 - n].net_make(net, nc, nd))
            na, nb = nc, nd

        na, nb = n3, n4

        # Draw components below centre
        for n in range(num_branches):

            if not (N & 1) and n == 0:
                sep = H[(N + 1) // 2]
            else:
                sep = H[(N + 1) // 2 + n] + H[(N + 1) // 2 - 1 + n]

            nc, nd =  net.node, net.node
            s.append('W %s %s; down=%s' % (na, nc, sep))
            s.append('W %s %s; down=%s' % (nb, nd, sep))
            s.append(self.args[(N + 1) // 2 + n].net_make(net, nc, nd))
            na, nb = nc, nd

        if n2 is None:
            n2 = net.node

        s.append('W %s %s; right=%s' % (n4, n2, self.wsep))
        return '\n'.join(s)


class Ser(ParSer):
    """Series class"""

    _operator = '+'

    def __init__(self, *args):

        _check_oneport_args(args)
        super(Ser, self).__init__()
        self.args = args

        for n, arg1 in enumerate(self.args):
            for arg2 in self.args[n + 1:]:
                if isinstance(arg1, I) and isinstance(arg2, I):
                    raise ValueError('Current sources connected in series'
                                     ' %s and %s' % (arg1, arg2))
                if isinstance(arg1, I) or isinstance(arg2, I):
                    # Should simplify by removing redundant component to
                    # save later grief with Norton or Thevenin transformation.
                    print('Warn: series connection with current source:'
                          ' %s and %s' % (arg1, arg2))

    @property
    def height(self):

        total = 0
        for arg in self.args:
            val = arg.height
            if val > total:
                total = val
        return total

    @property
    def width(self):

        total = 0
        for arg in self.args:
            total += arg.width
        return total + (len(self.args) - 1) * self.wsep

    def net_make(self, net, n1=None, n2=None):

        s = []
        if n1 is None:
            n1 = net.node
        for arg in self.args[:-1]:
            n3 = net.node
            s.append(arg.net_make(net, n1, n3))
            n1 = net.node
            s.append('W %s %s; right=%s' % (n3, n1, self.wsep))

        if n2 is None:
            n2 = net.node
        s.append(self.args[-1].net_make(net, n1, n2))
        return '\n'.join(s)


class R(OnePort):
    """Resistor"""

    def __init__(self, Rval):

        self.args = (Rval, )
        self.R = cExpr(Rval)
        self._Z = Zs.R(self.R)


class G(OnePort):
    """Conductance"""

    def __init__(self, Gval):

        self.args = (Gval, )
        self.G = cExpr(Gval)
        self._Z = 1 / Ys.G(self.G)

    def net_make(self, net, n1=None, n2=None):

        if n1 == None:
            n1 = net.node
        if n2 == None:
            n2 = net.node
        return 'R %s %s {%s}; right' % (n1, n2, 1 / self.G)


class L(OnePort):
    """Inductor

    Inductance Lval, initial current i0"""

    def __init__(self, Lval, i0=None):

        self.hasic = i0 is not None
        if i0 is None:
            i0 = 0

        if self.hasic:
            self.args = (Lval, i0)
        else:
            self.args = (Lval, )

        Lval = cExpr(Lval)
        i0 = cExpr(i0)
        self.L = Lval
        self.i0 = i0
        self._Z = Zs.L(Lval)
        self._Voc = Vsuper(-Vs(i0 * Lval))
        self.zeroic = self.i0 == 0 


class C(OnePort):
    """Capacitor

    Capacitance Cval, initial voltage v0"""

    def __init__(self, Cval, v0=None):

        self.hasic = v0 is not None
        if v0 is None:
            v0 = 0

        if self.hasic:
            self.args = (Cval, v0)
        else:
            self.args = (Cval, )

        Cval = cExpr(Cval)
        v0 = cExpr(v0)
        self.C = Cval
        self.v0 = v0
        self._Z = Zs.C(Cval)
        self._Voc = Vsuper(Vs(v0).integrate())
        self.zeroic = self.v0 == 0


class Y(OnePort):
    """General admittance."""

    def __init__(self, Yval):

        self.args = (Yval, )
        Yval = Ys(Yval)
        self._Z = 1 / Yval


class Z(OnePort):
    """General impedance."""

    def __init__(self, Zval):

        self.args = (Zval, )
        Zval = Zs(Zval)
        self._Z = Zval


class VoltageSource(OnePort):

    voltage_source = True
    netname = 'V'
    is_noisy = False


class sV(VoltageSource):
    """Arbitrary s-domain voltage source"""

    netkeyword = 's'

    def __init__(self, Vval):

        self.args = (Vval, )
        Vval = sExpr(Vval)
        self._Voc = Vsuper(Vs(Vval))


class V(VoltageSource):
    """Voltage source. If the expression contains s treat as s-domain
    voltage otherwise time domain.  A constant V is considered DC
    with an s-domain voltage V / s."""

    def __init__(self, Vval):

        self.args = (Vval, )
        self._Voc = Vsuper(Vval)

        
class Vstep(VoltageSource):
    """Step voltage source (s domain voltage of v / s)."""

    netkeyword = 'step'

    def __init__(self, v):

        self.args = (v, )
        v = cExpr(v)
        self._Voc = Vsuper(Vs(tExpr(v).laplace(), causal=True))
        self.v0 = v


class Vdc(VoltageSource):
    """DC voltage source (note a DC voltage source of voltage V has
    an s domain voltage of V / s)."""

    netkeyword = 'dc'
    
    def __init__(self, v):

        self.args = (v, )
        v = cExpr(v)
        self._Voc = Vsuper(Vconst(v, dc=True))
        self.v0 = v

    @property
    def voc(self):
        return self.v0


class Vac(VoltageSource):
    """AC voltage source."""

    netkeyword = 'ac'

    def __init__(self, V, phi=0):

        self.args = (V, phi)
        V = Expr(V)
        phi = Expr(phi)

        # Note, cos(-pi / 2) is not quite zero.

        self.omega = symbol('omega', real=True)
        self.v0 = V
        self.phi = phi
        self._Voc = Vsuper(Vphasor(self.v0 * exp(j * self.phi), ac=True,
                                  omega=self.omega))


    @property
    def voc(self):
        return self.v0 * cos(self.omega * t + self.phi)


class Vnoise(VoltageSource):
    """Noise voltage source."""

    netkeyword = 'noise'
    is_noisy = True

    def __init__(self, V):

        self.args = (V, )
        self._Voc = Vsuper(Vn(V))

        
class v(VoltageSource):
    """Arbitrary t-domain voltage source"""

    def __init__(self, vval):

        self.args = (vval, )
        Vval = tExpr(vval)
        self._Voc = Vsuper(Vval)


class CurrentSource(OnePort):

    current_source = True
    netname = 'I'
    is_noisy = False    

    
class sI(CurrentSource):
    """Arbitrary s-domain current source"""

    netkeyword = 's'

    def __init__(self, Ival):

        self.args = (Ival, )
        Ival = sExpr(Ival)
        self._Isc = Isuper(Is(Ival))


class I(CurrentSource):
    """Current source. If the expression contains s treat as s-domain
    current otherwise time domain.  A constant I is considered DC with
    an s-domain current I / s.

    """

    def __init__(self, Ival):

        self.args = (Ival, )
        self._Isc = Isuper(Ival)

            
class Istep(CurrentSource):
    """Step current source (s domain current of i / s)."""

    netkeyword = 'step'

    def __init__(self, i):

        self.args = (i, )
        i = cExpr(i)
        self._Isc = Isuper(Is(tExpr(i).laplace(), causal=True))        
        self.i0 = i


class Idc(CurrentSource):
    """DC current source (note a DC current source of current i has
    an s domain current of i / s)."""

    netkeyword = 'dc'
    
    def __init__(self, i):

        self.args = (i, )
        i = cExpr(i)
        self._Isc = Isuper(Iconst(i, dc=True))
        self.i0 = i

    @property
    def isc(self):
        return self.i0


class Iac(CurrentSource):
    """AC current source."""

    netkeyword = 'ac'

    def __init__(self, I, phi=0):

        self.args = (I, phi)
        I = Expr(I)
        phi = Expr(phi)

        self.omega = symbol('omega', real=True)
        self.i0 = I
        self.phi = phi
        self._Isc = Isuper(Iphasor(self.i0 * exp(j * self.phi), ac=True,
                                  omega=self.omega))

    @property
    def isc(self):
        return self.i0 * cos(self.omega * t + self.phi)


class Inoise(CurrentSource):
    """Noise current source."""

    netkeyword = 'noise'
    is_noisy = True

    def __init__(self, I):

        self.args = (I, )
        self._Isc = Isuper(In(I))

        
class i(CurrentSource):
    """Arbitrary t-domain current source"""

    def __init__(self, ival):

        self.args = (ival, )
        Ival = tExpr(ival)
        self._Isc = Isuper(Ival)


class Xtal(OnePort):
    """Crystal

    This is modelled as a series R, L, C circuit in parallel
    with C0 (a Butterworth van Dyke model).  Note,
    harmonic resonances are not modelled.
    """

    def __init__(self, C0, R1, L1, C1):

        self.C0 = cExpr(C0)
        self.R1 = cExpr(R1)
        self.L1 = cExpr(L1)
        self.C1 = cExpr(C1)

        N = self.expand()
        super(Xtal, self).__init__(N.Z, N.V)
        self.args = (C0, R1, L1, C1)

    def expand(self):

        return (R(self.R1) + L(self.L1) + C(self.C1)) | C(self.C0)


class FerriteBead(OnePort):
    """Ferrite bead (lossy inductor)

    This is modelled as a series resistor (Rs) connected
    to a parallel R, L, C network (Rp, Lp, Cp).
    """

    def __init__(self, Rs, Rp, Cp, Lp):

        self.Rs = cExpr(Rs)
        self.Rp = cExpr(Rp)
        self.Cp = cExpr(Cp)
        self.Lp = cExpr(Lp)

        N = self.expand()
        super(Xtal, self).__init__(N.Z, N.V)
        self.args = (Rs, Rp, Cp, Lp)

    def expand(self):

        return R(self.Rs) + (R(self.Rp) + L(self.Lp) + C(self.Cp))


# Import this at end to circumvent circular dependencies
from lcapy.twoport import Ladder, LSection, TSection
