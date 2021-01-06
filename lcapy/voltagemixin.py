from .quantity import Quantity

class VoltageMixin(Quantity):

    quantity = 'voltage'
    quantity_label = 'Voltage'
    units = 'V'
    is_voltage = True

    def cpt(self):
        from .oneport import V
        return V(self)

    
    
