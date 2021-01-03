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


class ExpressionClasses(dict):

    # These could be lazily created as required.  Each domain would
    # need to register itself first.
    
    def make(self, domainclass, quantities=None, suffix='Expression'):

        if quantities is None:
            quantities = list(quantityclasses.keys())
    
        classdict = {}
        classdict['undefined'] = domainclass

        for quantity in quantities:
            quantityclass = quantityclasses[quantity]
            quantityunits = quantityclass.units
            
            domainunits = domainclass.domain_units
            
            if domainunits == '':
                units = quantityunits
            else:
                units = '%s/%s' % (quantityunits, domainunits)            
            
            docstring = '%s-domain %s (units %s).' % (domainclass.domain_label,
                                                      quantity, units)
        
            name = domainclass.__name__.replace(suffix, quantityclass.quantity.capitalize())
            
            newclass = type(name, (quantityclass, domainclass),
                            {'__doc__': docstring,
                             'units': units})
            classdict[quantity] = newclass

        return classdict

    def add(self, domain, classes):
        self[domain] = classes

    def get_quantity(self, domain, quantity):
        """Return appropriate expression class for the specified domain and quantity."""

        if domain not in self:
            raise ValueError('Unknown domain %s' % domain)

        classes = self[domain]
        
        if quantity not in classes:
            raise ValueError('Unknown quantity %s for domain %s' % (quantity, domain))
        return classes[quantity]


expressionclasses = ExpressionClasses()
