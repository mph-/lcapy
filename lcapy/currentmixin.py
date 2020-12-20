class CurrentMixin(object):

    wrapper = 'current'
    is_current = True

    def cpt(self):
        from .oneport import I
        
        return I(self)

    
    
