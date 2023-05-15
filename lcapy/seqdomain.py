class SeqDomain(object):

    is_sequence = True

    def _class_get(self, quantity, domain):

        from .seqmap import seqmap
        return seqmap(quantity, domain)

    def _class_by_quantity(self, quantity, domain=None):

        if domain is None:
            domain = self.domain
        return self._class_get(quantity, domain)

    def _class_by_domain(self, domain):

        return self._class_get(self.quantity, domain)

    def as_quantity(self, quantity):

        if quantity == 'voltage':
            return self.as_voltage()
        elif quantity == 'current':
            return self.as_current()
        elif quantity == 'impedance':
            return self.as_impedance()
        elif quantity == 'admittance':
            return self.as_admittance()
        elif quantity == 'transfer':
            return self.as_transfer()
        elif quantity == 'power':
            return self.as_power()
        elif quantity == 'undefined':
            return self.as_expr()
        raise ValueError('Unknown quantity %s for %s' % (quantity, self))

    def as_domain(self, domain):

        if domain == 'discrete time':
            return self.as_discrete_time()
        elif domain == 'discrete fourier':
            return self.as_discrete_fourier()
        elif domain == 'Z':
            return self.as_Z()
        raise ValueError('Unknown domain %s for %s' % (domain, self))

    def as_voltage(self):
        return self._class_by_quantity('voltage')(self)

    def as_current(self):
        return self._class_by_quantity('current')(self)

    def as_admittance(self):
        return self._class_by_quantity('admittance')(self)

    def as_impedance(self):
        return self._class_by_quantity('impedance')(self)

    def as_transfer(self):
        return self._class_by_quantity('transfer')(self)

    def as_power(self):
        return self._class_by_quantity('power')(self)

    def change(self, arg, domain=None, units_scale=None, **assumptions):
        """Change sequence class."""

        if domain is None:
            domain = self.domain

        quantity = self.quantity

        cls = self._class_by_quantity(quantity, domain)
        ret = cls(arg, **assumptions)
        return ret
