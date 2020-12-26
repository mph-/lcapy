class VoltageMixin(object):

    quantity = 'voltage'
    quantity_label = 'Voltage'
    units = 'V'
    is_voltage = True

    def cpt(self):
        from .oneport import V
        return V(self)

    
    
