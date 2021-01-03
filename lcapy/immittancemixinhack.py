from .omegaexpr import AngularFourierDomainImpedance, AngularFourierDomainAdmittance
from .sexpr import LaplaceDomainImpedance, LaplaceDomainAdmittance
from .symbols import j, omega, jomega, s
from .immittancemixin import ImmittanceMixin

class ImmittanceMixinHack(ImmittanceMixin):

    @property
    def Yw(self):
        """Admittance  Y(omega)."""
        return AngularFourierDomainAdmittance(self.admittance._selectexpr(omega))

    @property
    def Zw(self):
        """Impedance  Z(omega)."""
        return AngularFourierDomainImpedance(self.impedance._selectexpr(omega))

    @property
    def Ys(self):
        """Generalized admittance  Y(s)."""
        return LaplaceDomainAdmittance(self.admittance._selectexpr(s))

    @property
    def Zs(self):
        """Generalized impedance  Z(s)."""
        return LaplaceDomainImpedance(self.impedance._selectexpr(s))
