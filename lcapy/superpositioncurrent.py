"""This module provides the SupoerPositionCurrent class.  This represents currents
as a superposition in different transform domains.

For example, the following expression is a superposition of a DC
component, an AC component, and a transient component:

I1 = SuperpositionCurrent('1 + 2 * cos(2 * pi * 3 * t) + 3 * u(t)')

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
from .symbols import j, omega
from .impedance import impedance
from .currentmixin import CurrentMixin


class SuperpositionCurrent(Superposition, CurrentMixin):

    def cpt(self):
        from .oneport import I
        # Perhaps should generate more specific components such as Idc?        
        return I(self.time())

    def _mul(self, x):
        if isinstance(x, (int, float)):
            return self.__scale__(x)

        if isinstance(x, Superposition):
            raise TypeError('Cannot multiply %s by %s. '
            'You need to extract a specific component, e.g., a.s * b.s' %
            (type(self).__name__, type(x).__name__))
        
        if not x.is_impedance:
            raise TypeError("Unsupported types for *: 'Current' and '%s'" %
                            type(x).__name__)
        if x.is_time_domain:
            raise TypeError("Cannot multiply by time-domain impedance.")
        
        obj = self
        if 't' in self and not x.is_constant_domain:
            obj = self.decompose()
        xs = x.laplace()
            
        new = SuperpositionVoltage()
        if 'dc' in obj:
            new += obj['dc'] * xs(0)
        for key in obj.ac_keys():
            new += obj[key] * xs(j * obj[key].omega)
        for key in obj.noise_keys():            
            new += obj[key] * xs(omega)
        if 's' in obj:
            new += obj['s'] * xs
        if 't' in obj:
            new += obj['t'] * x
        return new

    def _div(self, x):
        if isinstance(x, (int, float)):
            return self.__scale__(1 / x)

        if isinstance(x, Superposition):
            raise TypeError("""
            Cannot divide %s by %s.  You need to extract a specific component, e.g., a.s / b.s.  If you want a transfer function use a(s) / b(s)""" % (type(self).__name__, type(x).__name__))

        if not x.is_admittance:
            raise TypeError("Cannot divide '%s' by '%s'; require admittance" %
                            (type(self).__name__, type(x).__name__))

        Z = 1 / x
        return self * Z

    def __mul__(self, x):
        if False:
            raise ValueError('Cannot multiply superposition, need to convert to specific domain')
        return self._mul(x)
    
    def __rmul__(self, x):
        return self.__mul__(x)
    
    def __truediv__(self, x):
        if False:
            raise ValueError('Cannot divide superposition, need to convert to specific domain')        
        
        return self._div(x)

    
from .superpositionvoltage import SuperpositionVoltage
