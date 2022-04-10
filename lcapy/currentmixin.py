from .quantity import Quantity


class CurrentMixin(Quantity):

    quantity = 'current'
    quantity_label = 'Current'
    quantity_units = 'A'
    is_current = True
    is_signal = True

    def cpt(self):
        from .oneport import I, Idc, Iac

        if self.is_dc:
            return Idc(self.expr)
        elif self.is_phasor_domain:
            return Iac(self, omega=self.omega)

        return I(self.expr)
