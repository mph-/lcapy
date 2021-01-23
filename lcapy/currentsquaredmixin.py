from .quantity import Quantity

class CurrentSquaredMixin(Quantity):

    quantity = 'currentsquared'
    quantity_label = 'Current^2'
    quantity_units = 'A^2'
    is_currentsquared = True
