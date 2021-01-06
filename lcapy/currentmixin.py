from .quantity import Quantity

class CurrentMixin(Quantity):

    quantity = 'current'
    quantity_label = 'Current'
    units = 'A'
    is_current = True

    def cpt(self):
        from .oneport import I
        
        return I(self)

    
    
