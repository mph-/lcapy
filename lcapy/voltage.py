"""This module provides the Voltage class.  This represents voltages
as a superposition in different transform domains.

For example, the following expression is a superposition of a DC
component, an AC component, and a transient component:

V1 = Voltage('1 + 2 * cos(2 * pi * 3 * t) + 3 * u(t)')

V1(t) returns the time domain expression
V1(s) returns the Laplace domain expression
V1(omega) returns the Fourier domain expression with angular frequency
V1(f) returns the Fourier domain expression with linear frequency

V1.dc returns the DC component
V1.ac returns a dictionary of the AC components, keyed by the frequency
V1.transient returns the time-domain transient component

V1.is_dc returns True if a pure DC signal
V1.is_ac returns True if a pure AC signal
V1.is_transient returns True if a pure transient signal
V1.has_dc returns True if has a DC signal
V1.has_ac returns True if has an AC signal
V1.has_transient returns True if has a transient signal

Copyright 2019--2020 Michael Hayes, UCECE

"""

from .super import Super
from .sym import omega0sym

class Voltage(Super):

    def __init__(self, *args, **kwargs):
        self.type_map = {cExpr: Vconst, sExpr : Vs, noiseomegaExpr: Vnoisy,
                         fExpr: Vf, noisefExpr: Vfnoisy, omegaExpr: Vomega,
                         Phasor: Vphasor, tExpr : Vt}
        self.decompose_domains = {'s': Vs, 'ac': Vphasor, 'dc':
                                  Vconst, 'n': Vnoisy, 't': Vt}
        self.time_class = Vt
        self.laplace_class = Vs    

        super (Voltage, self).__init__(*args, **kwargs)
        
    def __rmul__(self, x):
        return self.__mul__(x)

    def __mul__(self, x):
        if isinstance(x, (int, float)):
            return self.__scale__(x)

        if isinstance(x, Super):
            raise TypeError('Cannot multiply %s by %s. '
            'You need to extract a specific component, e.g., a.s * b.s' %
            (type(self).__name__, type(x).__name__))

        if not isinstance(x, Admittance):
            raise TypeError("Unsupported types for *: 'Voltage' and '%s'" %
                            type(x).__name__)
        obj = self
        if x.has(s):
            obj = self.decompose()
        
        new = Current()
        if 'dc' in obj:
            # TODO, fix types
            new += Iconst(obj['dc'] * cExpr(x.jomega(0)))
        for key in obj.ac_keys():
            new += obj[key] * x.jomega(obj[key].omega)
        for key in obj.noise_keys():            
            new += obj[key] * x.jomega
        if 's' in obj:
            new += obj['s'] * x
        if 't' in obj:
            new += self['t'] * tExpr(x)
            
        return new

    def __div__(self, x):
        if isinstance(x, (int, float)):
            return self.__scale__(1 / x)

        if isinstance(x, Super):
            raise TypeError("""
            Cannot divide %s by %s.  You need to extract a specific component, e.g., a.s / b.s.  If you want a transfer function use a(s) / b(s)""" % (type(self).__name__, type(x).__name__))

        if not isinstance(x, Impedance):
            raise TypeError("Cannot divide '%s' by '%s'; require Impedance" %
                            (type(self).__name__, type(x).__name__))        

        return self * Admittance(1 / x)

    def __truediv__(self, x):
        return self.__div__(x)

    def cpt(self):
        from .oneport import V
        # Perhaps should generate more specific components such as Vdc?
        return V(self.time())

    
def Vname(name, kind, cache=False):

    if kind in ('s', 'laplace'):    
        return Vs(name + '(s)')
    elif kind in ('t', 'time'):
        return Vt(name.lower() + '(t)')
    elif kind in (omega0sym, omega0, 'ac'):
        return Vphasor(name + '(omega_0)')
    # Not caching is a hack to avoid conflicts of Vn1 with Vn1(s) etc.
    # when using subnetlists.  The alternative is a proper context
    # switch.  This would require every method to set the context.
    return expr(name, cache=cache)            


def Vtype(kind):
    
    if isinstance(kind, str) and kind[0] == 'n':
        return Vnoisy
    try:
        return {'ivp' : Vs, 's' : Vs, 'n' : Vnoisy,
                'ac' : Vphasor, 'dc' : Vconst, 't' : Vt, 'time' : Vt}[kind]
    except KeyError:
        return Vphasor


from .expr import expr    
from .cexpr import Vconst, Iconst, cExpr
from .fexpr import fExpr, Vf
from .omegaexpr import omegaExpr, Vomega
from .sexpr import sExpr, Vs
from .texpr import tExpr, Vt
from .noiseomegaexpr import noiseomegaExpr, Vnoisy
from .noisefexpr import noisefExpr, Vfnoisy
from .phasor import Vphasor, Phasor
from .impedance import Impedance
from .admittance import Admittance
from .omegaexpr import omegaExpr
from .symbols import s, omega0
from .current import Current
