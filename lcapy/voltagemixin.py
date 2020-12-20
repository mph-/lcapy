class VoltageMixin(object):

    wrapper = 'voltage'
    is_voltage = True

    def cpt(self):
        from .oneport import V
        return V(self)

    
    
