class ExprDomain(object):

    is_sequence = False

    def _class_get(self, quantity, domain):

        from .exprmap import exprmap
        return exprmap(quantity, domain)

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

        if domain == 'time':
            return self.as_time()
        elif domain == 'laplace':
            return self.as_laplace()
        elif domain == 'fourier':
            return self.as_fourier()
        elif domain == 'phasor':
            return self.as_phasor()
        elif domain == 'angular fourier':
            return self.as_angular_fourier()
        elif domain == 'frequency response':
            return self.as_frequency_response()
        elif domain == 'angular frequency response':
            return self.as_angular_frequency_response()
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

    def as_expr(self):
        return self

    def as_constant(self):
        from .cexpr import cexpr

        if not self.is_unchanging:
            raise ValueError('Expression %s is not constant' % self)
        return self._class_by_quantity(self.quantity)(self)(cexpr(self))

    def as_superposition(self):
        from .superpositionvoltage import SuperpositionVoltage
        from .superpositioncurrent import SuperpositionCurrent
        from .super import Superposition

        if self.is_voltage:
            return SuperpositionVoltage(self)
        elif self.is_current:
            return SuperpositionCurrent(self)
        elif self.is_undefined:
            return Superposition(self)
        raise ValueError(
            'Can only convert voltage or current to superposition')

    def change(self, arg, domain=None, units_scale=None, **assumptions):
        """Change expression class."""

        if domain is None:
            domain = self.domain

        if domain == 'constant':
            from .expr import expr

            # Allow changing of constants, e.g., V1 to 5 * t
            domain = expr(arg).domain

        quantity = self.quantity

        cls = self._class_by_quantity(quantity, domain)
        ret = cls(arg, **assumptions)

        if units_scale is not None:
            ret.units = self.units * units_scale
        return ret
