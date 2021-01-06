class Quantity(object):

    quantity = 'undefined'    
    is_voltage = False
    is_current = False
    is_impedance = False
    is_admittance = False
    is_transfer = False
    is_immitance = False
    is_always_causal = False


class UndefinedQuantity(Quantity):
    pass
