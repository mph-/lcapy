"""This module provides the Current class.  This represents currents
as a superposition in different transform domains.

For example, the following expression is a superposition of a DC
component, an AC component, and a transient component:

I1 = Current('1 + 2 * cos(2 * pi * 3 * t) + 3 * u(t)')

I1(t) returns the time domain expression
I1(s) returns the Laplace domain expression
I1(omega) returns the Fourier domain expression with angular frequency
I1(f) returns the Fourier domain expression with linear frequency

I1.dc returns the DC component
I1.ac returns a dictionary of the AC components, keyed by the frequency
I1.transient returns the time-domain transient component

I1.is_dc returns True if a pure DC signal
I1.is_ac returns True if a pure AC signal
I1.is_transient returns True if a pure transient signal
I1.has_dc returns True if has a DC signal
I1.has_ac returns True if has an AC signal
I1.has_transient returns True if has a transient signal

Copyright 2019--2020 Michael Hayes, UCECE

"""

from .super import Superposition
from .sym import omega0sym

class Current(Superposition):

    def __init__(self, *args, **kwargs):    

        self.type_map = {ConstantExpression: ConstantCurrent, LaplaceDomainExpression : LaplaceDomainCurrent, AngularFourierDomainNoiseExpression: AngularFourierDomainNoiseCurrent,
                         FourierDomainExpression: FourierDomainCurrent, FourierDomainNoiseExpression: FourierDomainNoiseCurrent, AngularFourierDomainExpression: AngularFourierDomainCurrent,
                         PhasorExpression: PhasorCurrent, TimeDomainExpression : TimeDomainCurrent}
        self.decompose_domains = {'s': LaplaceDomainCurrent, 'ac': PhasorCurrent, 'dc':
                                  ConstantCurrent, 'n': AngularFourierDomainNoiseCurrent, 't': TimeDomainCurrent}
        self.time_class = TimeDomainCurrent
        self.laplace_class = LaplaceDomainCurrent

        super (Current, self).__init__(*args, **kwargs)

    def __rmul__(self, x):
        return self.__mul__(x)
    
    def __mul__(self, x):
        if isinstance(x, (int, float)):
            return self.__scale__(x)

        if isinstance(x, Superposition):
            raise TypeError('Cannot multiply %s by %s. '
            'You need to extract a specific component, e.g., a.s * b.s' %
            (type(self).__name__, type(x).__name__))
        
        if not isinstance(x, Impedance):
            raise TypeError("Unsupported types for *: 'Current' and '%s'" %
                            type(x).__name__)
        obj = self
        if x.has(s):
            obj = self.decompose()

        new = Voltage()
        if 'dc' in obj:
            # TODO, fix types
            new += ConstantVoltage(obj['dc'] * ConstantExpression(x.jomega(0)))
        for key in obj.ac_keys():
            new += obj[key] * x.jomega(obj[key].omega)
        for key in obj.noise_keys():            
            new += obj[key] * x.jomega            
        if 's' in obj:
            new += obj['s'] * x
        if 't' in obj:
            new += obj['t'] * TimeDomainExpression(x)                        
        return new

    def __div__(self, x):
        if isinstance(x, (int, float)):
            return self.__scale__(1 / x)

        if isinstance(x, Superposition):
            raise TypeError("""
            Cannot divide %s by %s.  You need to extract a specific component, e.g., a.s / b.s.  If you want a transfer function use a(s) / b(s)""" % (type(self).__name__, type(x).__name__))

        if not isinstance(x, Admittance):
            raise TypeError("Cannot divide '%s' by '%s'; require Admittance" %
                            (type(self).__name__, type(x).__name__))

        return self * Impedance(1 / x)

    def __truediv__(self, x):
        return self.__div__(x)

    def cpt(self):
        from .oneport import I
        # Perhaps should generate more specific components such as Idc?        
        return I(self.time())

    
def Iname(name, kind, cache=False):
    
    if kind in ('s', 'laplace'):
        return LaplaceDomainCurrent(name + '(s)')
    elif kind in ('t', 'time'):
        return TimeDomainCurrent(name.lower() + '(t)')
    elif kind in (omega0sym, omega0, 'ac'):    
        return PhasorCurrent(name + '(omega_0)')
    return expr(name, cache=cache)            


def Itype(kind):
    if isinstance(kind, str) and kind[0] == 'n':
        return AngularFourierDomainNoiseCurrent
    try:
        return {'ivp' : LaplaceDomainCurrent, 's' : LaplaceDomainCurrent, 'n' : AngularFourierDomainNoiseCurrent,
                'ac' : PhasorCurrent, 'dc' : ConstantCurrent, 't' : TimeDomainCurrent, 'time' : TimeDomainCurrent}[kind]
    except KeyError:
        return PhasorCurrent
    

from .expr import expr
from .cexpr import ConstantVoltage, ConstantCurrent, ConstantExpression        
from .fexpr import FourierDomainExpression, FourierDomainCurrent
from .omegaexpr import AngularFourierDomainExpression, AngularFourierDomainCurrent
from .sexpr import LaplaceDomainExpression, LaplaceDomainCurrent
from .texpr import TimeDomainExpression, TimeDomainCurrent
from .noiseomegaexpr import AngularFourierDomainNoiseExpression, AngularFourierDomainNoiseCurrent
from .noisefexpr import FourierDomainNoiseExpression, FourierDomainNoiseCurrent
from .phasor import PhasorCurrent, PhasorExpression
from .impedance import Impedance
from .admittance import Admittance
from .omegaexpr import AngularFourierDomainExpression
from .symbols import s, omega0
from .voltage import Voltage

