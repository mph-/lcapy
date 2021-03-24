"""This module provides the ExpressionClasses class, for dynamically
creating expression classes, and a factory to select the appropriate
class for a domain and quantity.

Copyright 2021 Michael Hayes, UCECE

"""

from .voltagemixin import VoltageMixin
from .currentmixin import CurrentMixin
from .admittancemixin import AdmittanceMixin
from .impedancemixin import ImpedanceMixin
from .transfermixin import TransferMixin
from .voltagesquaredmixin import VoltageSquaredMixin
from .currentsquaredmixin import CurrentSquaredMixin
from .admittancesquaredmixin import AdmittanceSquaredMixin
from .impedancesquaredmixin import ImpedanceSquaredMixin
from .powermixin import PowerMixin
from .units import units, u as uu
from sympy import sqrt
from sympy import S


quantityclasses = {'voltage': VoltageMixin,
                   'current': CurrentMixin,
                   'admittance': AdmittanceMixin,
                   'impedance': ImpedanceMixin,
                   'transfer': TransferMixin,
                   'voltagesquared': VoltageSquaredMixin,
                   'currentsquared': CurrentSquaredMixin,
                   'admittancesquared': AdmittanceSquaredMixin,
                   'impedancesquared': ImpedanceSquaredMixin,
                   'power': PowerMixin}

units_mapping = {
    '': S.One,
    'V': uu.volt, 'A': uu.ampere, 
    'V/Hz': uu.volt / uu.Hz, 'A/Hz': uu.ampere / uu.Hz,
    'V/sqrt(Hz)': uu.volt / sqrt(uu.Hz), 'A/sqrt(Hz)': uu.ampere / sqrt(uu.Hz),
    'ohm': uu.ohm, 'S': uu.siemens,
    'ohm/s': uu.ohm / uu.s, 'S/s': uu.siemens / uu.s,
    'ohm^2/s^2': (uu.ohm / uu.s)**2, 'S^2/s^2': (uu.siemens / uu.s)**2,    
    'V^2': uu.volt**2, 'A^2': uu.ampere**2, 
    'V^2/Hz^2': (uu.volt / uu.Hz)**2, 'A^2/Hz^2': (uu.ampere / uu.Hz)**2,
    'ohm^2': uu.ohm**2, 'S^2': uu.siemens**2,
    'W': uu.watt, '/s': 1 / uu.s}


class ExpressionClassBuilder(dict):

    def __init__(self, domain, domainclass1, domainclass2=None, quantities=None):

        self.domain = domain
        self.domainclass1 = domainclass1
        self.domainclass2 = domainclass2
        self.quantities = quantities

    def make1(self, quantity, domainclass):

        if quantity == 'undefined':
            self[quantity] = domainclass
            return domainclass

        quantityclass = quantityclasses[quantity]
        quantityunits = quantityclass.quantity_units
            
        unitsstring = quantityunits            
        if quantity in ('voltage', 'current'):
            if (domainclass.is_laplace_domain or
                domainclass.is_fourier_domain or
                domainclass.is_angular_fourier_domain):
                unitsstring = '%s/Hz' % quantityunits
        elif quantity in ('voltagesquared', 'currentsquared'):
            if (domainclass.is_laplace_domain or
                domainclass.is_fourier_domain or
                domainclass.is_angular_fourier_domain):                
                unitsstring = '%s/Hz^2' % quantityunits
        elif quantity in ('impedance', 'admittance', 'transfer'):
            if domainclass.is_time_domain:
                unitsstring = '%s/s' % quantityunits
        elif quantity in ('admittancesquared', 'impedancesquared'):
            if domainclass.is_time_domain:
                unitsstring = '%s/s^2' % quantityunits

        # TODO: perhaps rename transfer function in time domain to
        # impulse response?
        docstring = '%s-domain %s (units %s).' % (domainclass.domain_label,
                                                  quantity, unitsstring)

        suffix = 'Expression'
        name = domainclass.__name__.replace(suffix, quantityclass.quantity.capitalize())
            
        newclass = type(name, (quantityclass, domainclass),
                            {'__doc__': docstring,
                             '_default_units': units_mapping[unitsstring]})
        self[quantity] = newclass

        #print('Created %s %s' % (self.domain, quantity))
        return newclass
        
    def make(self, quantity):
        
        if self.quantities is None:
            return self.make1(quantity, self.domainclass1)

        if quantity in self.quantities:
            return self.make1(quantity, self.domainclass1)

        return self.make1(quantity, self.domainclass2)        

    def __getitem__(self, quantity):

        if quantity in self:
            return super(ExpressionClassBuilder, self).__getitem__(quantity)

        return self.make(quantity)
    

class ExpressionClasses(dict):                 
    
    def register(self, domain, domainclass1, domainclass2=None, quantities=None):

        self[domain] = ExpressionClassBuilder(domain, domainclass1, domainclass2, quantities)
        return self[domain]
                 
    def get_quantity(self, domain, quantity):

        return self[domain][quantity]


expressionclasses = ExpressionClasses()    
