from .quantity import Quantity

class PhaseMixin(Quantity):

    quantity = 'phase'
    quantity_label = 'Phase'
    units = 'rad'
    is_phase = True
