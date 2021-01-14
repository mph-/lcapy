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
from .phasemixin import PhaseMixin


quantityclasses = {'voltage': VoltageMixin,
                   'current': CurrentMixin,
                   'admittance': AdmittanceMixin,
                   'impedance': ImpedanceMixin,
                   'transfer': TransferMixin,
                   'voltagesquared': VoltageSquaredMixin,
                   'currentsquared': CurrentSquaredMixin,
                   'admittancesquared': AdmittanceSquaredMixin,
                   'impedancesquared': ImpedanceSquaredMixin,
                   'power': PowerMixin,
                   'phase': PhaseMixin}


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

            units = quantityunits            
            if quantity in ('voltage', 'current'):
                if not (domainclass.is_time_domain or
                        domainclass.is_discrete_time_domain or
                        domainclass.is_constant_domain or
                        domainclass.is_phasor_domain):                
                    units = '%s/%s' % (quantityunits, domainunits)
            elif quantity in ('voltagesquared', 'currentsquared'):
                if not (domainclass.is_time_domain or
                        domainclass.is_discrete_time_domain or
                        domainclass.is_constant_domain or
                        domainclass.is_phasor_domain):                
                    units = '%s/%s^2' % (quantityunits, domainunits)
            elif quantity in ('impedance', 'admittance'):
                if domainclass.is_time_domain or domainclass.is_discrete_time_domain:
                    units = '%s/%s' % (quantityunits, domainunits)
            elif quantity in ('admittancesquared', 'impedancesquared'):
                if domainclass.is_time_domain or domainclass.is_discrete_time_domain:
                    units = '%s/%s^2' % (quantityunits, domainunits)

            # FIXME:  The units of squared quantities are incorrect under transformation
            # to another domain.  For example, v1(t) * v2(t) has units V^2,
            # V1(f) * V2(f) has units (V/Hz)^2, but (v1(t) * v2(t))(f) has units V^2/Hz,
            # and (V1(f) * V2(f))(t) has units V^2/Hz.
                        
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
