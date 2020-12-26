class CurrentMixin(object):

    quantity = 'current'
    quantity_label = 'Current'
    units = 'A'
    is_current = True

    def cpt(self):
        from .oneport import I
        
        return I(self)

    
    
