from .quantity import Quantity

class VoltageMixin(Quantity):

    quantity = 'voltage'
    quantity_label = 'Voltage'
    units = 'V'
    is_voltage = True

    def cpt(self):
        from .oneport import V, Vac, Vdc

        if self.is_dc:
            return Vdc(self.expr)
        elif self.is_ac:
            p = self.as_phasor()
            return Vac(p.expr, omega=p.omega)
        
        return V(self.expr)

    
    
