from .quantity import Quantity

class VoltageSquaredMixin(Quantity):

    quantity = 'voltagesquared'
    quantity_label = 'Voltage^2'
    units = 'V^2'
    is_voltagesquared = True


    
