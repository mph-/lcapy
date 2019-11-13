from .cexpr import cExpr
from .sexpr import sExpr
from .omegaexpr import omegaExpr
from .symbols import j, omega, jomega, s
from .sym import omegasym


class ZY(sExpr):
    
    def __init__(self, val, kind=None, causal=True, **assumptions):

        super(ZY, self).__init__(val, causal=causal, **assumptions)
        self.kind = kind

    @property
    def pexpr(self):

        kind = self.kind
        if kind is None:
            # Default to printing impedance as Z(omega)
            kind = omega
        
        return self.select(kind).expr
        
    @property            
    def jomega(self):
        return self(jomega)

    @property    
    def s(self):
        return self(s)    

    def select(self, kind=None):

        if kind is None:
            kind = self.kind

        if kind in ('s', 'ivp', 'super', 'laplace'):
            return self
        elif kind in ('dc', 'time'):
            return cExpr(self.subs(0))
        elif isinstance(kind, str) and kind[0] == 'n':
            return self(jomega)
        elif kind in (omegasym, omega, 'ac'):
            return self(jomega)
        return omegaExpr(self.subs(j * kind))    

    @property
    def R(self):
        """Resistance."""
        return self.Zw.real

    @property
    def X(self):
        """Reactance."""
        return self.Zw.imag

    @property
    def G(self):
        """Conductance."""
        return self.Yw.real

    @property
    def B(self):
        """Susceptance."""
        return -self.Yw.imag    

class ZZ(ZY):
    """Generic impedance class."""

    @property
    def Yw(self):
        return 1 / self.Zw

    @property
    def Zw(self):
        return self.jomega

    def cpt(self):
        from .oneport import R, C, L, Z

        if self.is_number or self.is_dc:
            return R(self.expr)

        z = self * s

        if z.is_number:
            return C((1 / z).expr)

        z = self / s

        if z.is_number:
            return L(z.expr)

        return Z(self)
    

class YY(ZY):    
    """Generic admittance class."""

    @property
    def Yw(self):
        return self.jomega

    @property
    def Zw(self):
        return 1 / self.Yw

    def cpt(self):
        from .oneport import G, C, L, Y

        if self.is_number or self.is_dc:
            return G(self.expr)

        y = self * s

        if y.is_number:
            return L((1 / y).expr)

        y = self / s

        if y.is_number:
            return C(y.expr)

        return Y(self)    

    
