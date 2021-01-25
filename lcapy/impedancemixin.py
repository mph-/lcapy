from .expr import expr
from .immittancemixin import ImmittanceMixin
from .quantity import Quantity

class ImpedanceMixin(Quantity, ImmittanceMixin):

    quantity = 'impedance'
    quantity_label = 'Impedance'
    quantity_units = 'ohm'
    is_impedance = True
    is_immittance = True
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

        from .admittance import admittance

        return admittance(1 / self)

    @property
    def Z(self):
        """Impedance."""
        return self
    
    def __rtruediv__(self, x):
        """Reverse true divide"""

        x = expr(x)
        if x.is_constant:
            from .admittance import admittance            
            return admittance(x.expr / self.expr)
        return super(ImpedanceMixin, self).__rtruediv__(x)
    
    def cpt(self):
        from .oneport import R, C, L, Z
        from .symbols import s        

        if self.is_number or self.is_dc:
            return R(self.expr)

        z = self.laplace() * s

        if z.is_number:
            return C((1 / z).expr)

        z = self.laplace() / s

        if z.is_number:
            return L(z.expr)

        return Z(self)

