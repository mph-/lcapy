"""This module provides the SuperpositionVoltage class.  This represents voltages
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

from .super import Superposition
from .sym import omega0sym

class SuperpositionVoltage(Superposition):

    def __init__(self, *args, **kwargs):
        self.type_map = {ConstantExpression: ConstantVoltage,
                         LaplaceDomainExpression : LaplaceDomainVoltage,
                         AngularFourierDomainNoiseExpression: AngularFourierDomainNoiseVoltage,
                         FourierDomainExpression: FourierDomainVoltage,
                         FourierDomainNoiseExpression: FourierDomainNoiseVoltage,
                         AngularFourierDomainExpression: AngularFourierDomainVoltage,
                         PhasorDomainExpression: PhasorDomainVoltage,
                         TimeDomainExpression : TimeDomainVoltage}
        self.decompose_domains = {'s': LaplaceDomainVoltage,
                                  'ac': PhasorDomainVoltage,
                                  'dc': ConstantVoltage,
                                  'n': AngularFourierDomainNoiseVoltage,
                                  't': TimeDomainVoltage}
        self.time_class = TimeDomainVoltage
        self.laplace_class = LaplaceDomainVoltage    

        super (SuperpositionVoltage, self).__init__(*args, **kwargs)
        
    def __rmul__(self, x):
        return self.__mul__(x)

    def __mul__(self, x):
        if isinstance(x, (int, float)):
            return self.__scale__(x)

        if isinstance(x, Superposition):
            raise TypeError('Cannot multiply %s by %s. '
            'You need to extract a specific component, e.g., a.s * b.s' %
            (type(self).__name__, type(x).__name__))

        if not x.is_admittance:
            raise TypeError("Unsupported types for *: 'Voltage' and '%s'" %
                            type(x).__name__)
        obj = self
        if x.has(s):
            obj = self.decompose()
        
        new = SuperpositionCurrent()
        if 'dc' in obj:
            # TODO, fix types
            new += ConstantCurrent(obj['dc'] * ConstantExpression(x.jomega(0)))
        for key in obj.ac_keys():
            new += obj[key] * x.jomega(obj[key].omega)
        for key in obj.noise_keys():            
            new += obj[key] * x.jomega
        if 's' in obj:
            new += obj['s'] * x
        if 't' in obj:
            new += self['t'] * TimeDomainExpression(x)
            
        return new

    def __div__(self, x):
        if isinstance(x, (int, float)):
            return self.__scale__(1 / x)

        if isinstance(x, Superposition):
            raise TypeError("""
            Cannot divide %s by %s.  You need to extract a specific component, e.g., a.s / b.s.  If you want a transfer function use a(s) / b(s)""" % (type(self).__name__, type(x).__name__))

        if not x.is_impedance:
            raise TypeError("Cannot divide '%s' by '%s'; require impedance" %
                            (type(self).__name__, type(x).__name__))        

        return self * admittance(1 / x)

    def __truediv__(self, x):
        return self.__div__(x)

    def cpt(self):
        from .oneport import V
        # Perhaps should generate more specific components such as Vdc?
        return V(self.time())

    
def Vname(name, kind, cache=False):

    if kind in ('s', 'laplace'):    
        return LaplaceDomainVoltage(name + '(s)')
    elif kind in ('t', 'time'):
        return TimeDomainVoltage(name.lower() + '(t)')
    elif kind in (omega0sym, omega0, 'ac'):
        return PhasorDomainVoltage(name + '(omega_0)')
    # Not caching is a hack to avoid conflicts of Vn1 with Vn1(s) etc.
    # when using subnetlists.  The alternative is a proper context
    # switch.  This would require every method to set the context.
    return expr(name, cache=cache)            


def Vtype(kind):
    
    if isinstance(kind, str) and kind[0] == 'n':
        return AngularFourierDomainNoiseVoltage
    try:
        return {'ivp' : LaplaceDomainVoltage, 's' : LaplaceDomainVoltage, 'n' : AngularFourierDomainNoiseVoltage,
                'ac' : PhasorDomainVoltage, 'dc' : ConstantVoltage, 't' : TimeDomainVoltage, 'time' : TimeDomainVoltage}[kind]
    except KeyError:
        return PhasorDomainVoltage


from .expr import expr    
from .cexpr import ConstantVoltage, ConstantCurrent, ConstantExpression
from .fexpr import FourierDomainExpression, FourierDomainVoltage
from .omegaexpr import AngularFourierDomainExpression, AngularFourierDomainVoltage
from .sexpr import LaplaceDomainExpression, LaplaceDomainVoltage
from .texpr import TimeDomainExpression, TimeDomainVoltage
from .noiseomegaexpr import AngularFourierDomainNoiseExpression, AngularFourierDomainNoiseVoltage
from .noisefexpr import FourierDomainNoiseExpression, FourierDomainNoiseVoltage
from .phasor import PhasorDomainVoltage, PhasorDomainExpression
from .impedance import impedance
from .admittance import admittance
from .omegaexpr import AngularFourierDomainExpression
from .symbols import s, omega0
from .superposition_current import SuperpositionCurrent
