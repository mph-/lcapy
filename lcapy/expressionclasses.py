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
        quantityunits = quantityclass.units
            
        domainunits = domainclass.domain_units

        units = quantityunits            
        if quantity in ('voltage', 'current'):
            if (domainclass.is_laplace_domain or
                domainclass.is_fourier_domain or
                domainclass.is_angular_fourier_domain):
                units = '%s/%s' % (quantityunits, domainunits)
        elif quantity in ('voltagesquared', 'currentsquared'):
            if (domainclass.is_laplace_domain or
                domainclass.is_fourier_domain or
                domainclass.is_angular_fourier_domain):                
                units = '%s/%s^2' % (quantityunits, domainunits)
        elif quantity in ('impedance', 'admittance'):
            if domainclass.is_time_domain:
                units = '%s/%s' % (quantityunits, domainunits)
        elif quantity in ('admittancesquared', 'impedancesquared'):
            if domainclass.is_time_domain:
                units = '%s/%s^2' % (quantityunits, domainunits)

        # FIXME:  The units of squared quantities are incorrect under transformation
        # to another domain.  For example, v1(t) * v2(t) has units V^2,
        # V1(f) * V2(f) has units (V/Hz)^2, but (v1(t) * v2(t))(f) has units V^2/Hz,
        # and (V1(f) * V2(f))(t) has units V^2/Hz.
                        
        docstring = '%s-domain %s (units %s).' % (domainclass.domain_label,
                                                  quantity, units)

        suffix = 'Expression'
        name = domainclass.__name__.replace(suffix, quantityclass.quantity.capitalize())
            
        newclass = type(name, (quantityclass, domainclass),
                            {'__doc__': docstring,
                             'units': units})
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
