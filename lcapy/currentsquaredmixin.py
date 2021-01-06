from .quantity import Quantity

class CurrentSquaredMixin(Quantity):

    quantity = 'currentsquared'
    quantity_label = 'Current^2'
    units = 'A^2'
    is_currentsquared = True
