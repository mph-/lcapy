"""This module provides the Units class for simplification of units.
It should be rolled into SymPy.  It can perform simplification of
units, e.g., volts / amperes -> ohms.

Copyright 2020--2021 Michael Hayes, UCECE

"""

import sympy.physics.units as u
from sympy.physics.units.systems.si import dimsys_SI
from sympy.physics.units.systems import SI 
from sympy.physics.units import UnitSystem, Quantity
from sympy import sqrt
from sympy import S


dB = Quantity('dB', 'dB')

units_mapping = {
    '': S.One,
    'V': u.volt, 'A': u.ampere, 
    'V/Hz': u.volt / u.Hz, 'A/Hz': u.ampere / u.Hz,
    'V/rad/s': u.volt / (u.rad / u.s) , 'A/rad/s': u.ampere / (u.rad / u.s),    
    'V/sqrt(Hz)': u.volt / sqrt(u.Hz), 'A/sqrt(Hz)': u.ampere / sqrt(u.Hz),
    'ohm': u.ohm, 'S': u.siemens,
    'ohm/s': u.ohm / u.s, 'S/s': u.siemens / u.s,
    'ohm^2/s^2': (u.ohm / u.s)**2, 'S^2/s^2': (u.siemens / u.s)**2,    
    'V^2': u.volt**2, 'A^2': u.ampere**2, 
    'V^2/Hz^2': (u.volt / u.Hz)**2, 'A^2/Hz^2': (u.ampere / u.Hz)**2,
    'V^2/(rad/s)^2': (u.volt / (u.rad / u.s)) ** 2, 'A^2/(rad/s^2)': (u.ampere / (u.rad / u.s))**2,        
    'ohm^2': u.ohm**2, 'S^2': u.siemens**2,
    'W': u.watt}


class Units(object):

    def __init__(self, unit_system="SI"):

        self.unit_system = UnitSystem.get_unit_system(unit_system)
        self.dim_sys = self.unit_system.get_dimension_system()
        self._mapping = {}
        
        for i in u.__dict__:          
            unit = getattr(u, i)
            if not isinstance(unit, u.Quantity):
                continue

            key = self._makekey(unit)            

            # Use earlier defined units
            
            if key not in self._mapping:
                self._mapping[key] = unit

        # Remove entry for no units.
        self._mapping.pop(self._makekey(1))

        # Add entry for S * ohm, etc.
        key = (None, ) * len(key)
        self._mapping[key] = S.One             

    def _get_dependencies(self, unit):
            
        dim = self.unit_system.get_dimensional_expr(unit)
        return self.dim_sys.get_dimensional_dependencies(dim)
            
    def _makekey(self, unit):
            
        deps = self._get_dependencies(unit)
        key = tuple([deps.get(str(dim.name)) for dim in self.dim_sys.base_dims])
        return key
            
    def simplify_units(self, unit):

        key = self._makekey(unit)
        if not key in self._mapping:
            return unit
        result = self._mapping[key]

        # V s or Wb?  In the context of Laplace transforms, V s makes more
        # sense since the Laplace domain voltage has units (V / rad / s).
        # However, for magnetic field strength, Wb makes more sense.  Since
        # this is for circuit analysis we plump for V s.
        if result.has(u.webers):
            result = result.replace(u.webers, u.volt * u.s)
        return result

    def simplify(self, expr):

        value, unit = self.as_value_unit(expr)
        return value * self.simplify_units(unit)

    def as_value_unit(self, expr):
        return as_value_unit(expr)

    def by_name(self, name):
        return units_mapping[name]

    
def as_value_unit(expr):

    if isinstance(expr, u.Quantity):
        return 1, expr

    if not expr.has(u.Quantity):
        return expr, 1

    if expr.is_Pow and expr.args[1] == -1:
        value, unit = as_value_unit(expr.args[0])
        return S.one / value, S.one / unit
    
    defs = {x: 1 for x in expr.args if not x.has(u.Quantity)}
    unit = expr.subs(defs)

    value = expr / unit
    if value.has(u.Quantity):
        # FIXME: This function only works for something like 42 * volt or 42 * amp * ohm.
        # It fails for 4 * amp * 2 * ohm + 42 * volt.
        raise ValueError('Expression not of form value * units: %s' % expr)
    
    return value, unit


units = Units()
