from .expr import expr
from .immittancemixin import ImmittanceMixin
from .quantity import Quantity

class AdmittanceMixin(Quantity, ImmittanceMixin):

    quantity = 'admittance'
    quantity_label = 'Admittance'
    quantity_units = 'S'
    is_admittance = True
    is_immitance = True
    is_ratio = True    

    # Immittances derived from a realisable circuit will be causal but
    # non-causal immittances can also be constructed.  So this might
    # disappear.  An example non-causal impulse response is z(t) = R delta(t + T)
    # with an impedance Z(s) = R exp(s * T).  This has
    # a real part R * exp(re(s) * T) * cos(T * im(s))
    # and imaginary part R * exp(re(s) * T) * sin(T * im(s))
    
    is_always_causal = True    

    @property
    def Y(self):
        """Admittance."""
        return self

    @property
    def Z(self):
        """Impedance."""

        from .impedance import impedance
        
        return impedance(1 / self)

    def __rtruediv__(self, x):
        """Reverse true divide"""

        x = expr(x)
        if x.is_constant:
            from .impedance import impedance            
            return impedance(x.expr / self.expr)            
        return super(AdmittanceMixin, self).__rtruediv__(x)    

    def cpt(self):
        from .oneport import G, C, L, Y
        from .symbols import s

        if self.is_number or self.is_dc:
            return G(self.expr)

        y = self.laplace() * s

        if y.is_number:
            return L((1 / y).expr)

        y = self.laplace() / s

        if y.is_number:
            return C(y.expr)

        return Y(self)
