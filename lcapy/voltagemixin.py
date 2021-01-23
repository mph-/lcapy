from .quantity import Quantity

class VoltageMixin(Quantity):

    quantity = 'voltage'
    quantity_label = 'Voltage'
    quantity_units = 'V'
    is_voltage = True

    def cpt(self):
        from .oneport import V, Vac, Vdc

        if self.is_dc:
            return Vdc(self.expr)
        elif self.is_phasor_domain:
            return Vac(self, omega=self.omega)
        
        return V(self.expr)

    
    
