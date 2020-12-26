from .expr import expr
from .immittancemixin import ImmittanceMixin

class AdmittanceMixin(ImmittanceMixin):

    quantity = 'admittance'
    quantity_label = 'Admittance'
    units = 'S'
    is_admittance = True

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
