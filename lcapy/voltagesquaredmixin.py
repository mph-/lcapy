from .quantity import Quantity

class VoltageSquaredMixin(Quantity):

    quantity = 'voltagesquared'
    quantity_label = 'Voltage^2'
    quantity_units = 'V^2'
    is_voltagesquared = True


    
