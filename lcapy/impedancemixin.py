from .expr import expr
from .immittancemixin import ImmittanceMixin

class ImpedanceMixin(ImmittanceMixin):

    wrapper = 'impedance'
    is_impedance = True

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

