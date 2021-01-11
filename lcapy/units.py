"""This module provides the Units class for simplification of units.
It should be rolled into SymPy.  It can perform simplification of
units, e.g., volts / amperes -> ohms.

Copyright 2020--2021 Michael Hayes, UCECE

"""

import sympy.physics.units as u
from sympy.physics.units.systems.si import dimsys_SI
from sympy.physics.units.systems import SI 
from sympy.physics.units import UnitSystem


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
        return self._mapping[key]

    def simplify(self, expr):

        value, unit = self.as_value_unit(expr)
        return value * self.simplify_units(unit)

    def as_value_unit(self, expr):
        return as_value_unit(expr)

    
def as_value_unit(expr):

    if isinstance(expr, u.Quantity):
        return 1, expr

    if not expr.has(u.Quantity):
        return expr, 1
    
    defs = {x: 1 for x in expr.args if not x.has(u.Quantity)}
    unit = expr.subs(defs)

    return expr / unit, unit


units = Units()
