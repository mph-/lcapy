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


quantities = {'voltage': VoltageMixin,
              'current': CurrentMixin,
              'admittance': AdmittanceMixin,
              'impedance': ImpedanceMixin,
              'transfer': TransferMixin,
              'voltagesquared': VoltageSquaredMixin,
              'currentsquared': CurrentSquaredMixin,
              'admittancesquared': AdmittanceSquaredMixin,
              'impedancesquared': ImpedanceSquaredMixin,
              'power': PowerMixin}
