class Quantity(object):

    quantity = ''    
    is_voltage = False
    is_current = False
    is_impedance = False
    is_admittance = False
    is_transfer = False
    is_immitance = False
    is_always_causal = False
    is_undefined = False


class UndefinedQuantity(Quantity):

    quantity = 'undefined'        
    is_undefined = True    

