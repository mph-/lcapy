"""This module provides the Quantity class, for tracking of an expression quantity.

Copyright 2020--2021 Michael Hayes, UCECE

"""

class Quantity(object):

    quantity = ''    
    is_power = False
    is_voltage = False
    is_current = False
    is_impedance = False
    is_admittance = False
    is_transfer = False
    is_immitance = False
    is_always_causal = False
    is_undefined = False
    is_voltagesquared = False
    is_currentsquared = False
    is_impedancesquared = False
    is_admittancesquared = False
    part = ''


class UndefinedQuantity(Quantity):

    quantity = 'undefined'        
    is_undefined = True    

