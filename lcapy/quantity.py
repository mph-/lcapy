"""This module provides the Quantity class, for tracking of an expression quantity.

Copyright 2020--2021 Michael Hayes, UCECE

"""


class Quantity(object):
    """This is the base class for VoltageMixin, etc.  An Lcapy quantity is
    not a true quantity but a collection of related quantities, for
    example, voltage (with units V) and voltage spectral density (with units V/Hz)."""

    quantity = ''
    is_power = False
    is_voltage = False
    is_current = False
    is_impedance = False
    is_admittance = False
    is_transfer = False
    is_immittance = False
    is_always_causal = False
    is_undefined = False
    is_voltagesquared = False
    is_currentsquared = False
    is_impedancesquared = False
    is_admittancesquared = False
    is_ratio = False
    is_signal = False
    is_squared = False
    part = ''


class UndefinedQuantity(Quantity):

    quantity = 'undefined'
    is_undefined = True
