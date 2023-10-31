"""This module provides the SubNetlist class.  It provides circuit
analysis for a single domain quantity.  It is for internal use only.
The Netlist class is for users.

Copyright 2020--2023 Michael Hayes, UCECE

"""

from .symbols import omega
from .state import state
from .mna import MNA
from .netlistmixin import NetlistMixin
from .netlistsimplifymixin import NetlistSimplifyMixin
from .netfile import NetfileMixin


class SubNetlist(NetlistMixin, NetlistSimplifyMixin, NetfileMixin):
    """This is a representation of a netlist for a particular
    transformation domain, such as ac, dc, transient, or noise.  It is
    for internal use only.  Unlike Netlist, SubNetlist is not mutable.
    """

    def __new__(cls, netlist, kind):

        kinds = ('t', 'dc', 's', 'time', 'ivp', 'laplace')
        if not isinstance(kind, str) or kind[0] == 'n':
            pass
        elif kind not in kinds:
            raise ValueError('Expected one of %s for kind, got %s' %
                             (', '.join(kinds), kind))

        obj = netlist.select(kind=kind)
        # Need own context to avoid conflicts with Vn1 and Vn1(s), etc.
        obj.context = state.new_context()
        obj.kind = kind
        obj.__class__ = cls
        obj._analysis = obj.analyse()
        obj.solver_method = netlist.solver_method
        return obj

    def __init__(self, netlist, kind):
        """ kind can be 't', 'dc', 's', 'time', 'ivp', 'n*' or omega,
        where 'n*' is a noise identifer and omega is an angular frequency."""

        self.mna = MNA(self, self.solver_method)

    def get_I(self, name):
        """Current through component"""

        return self.mna.Idict[name].canonical()

    def get_i(self, name):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def get_Vd(self, Np, Nm=None):
        """Voltage drop between nodes"""

        return (self.mna.Vdict[Np] - self.mna.Vdict[Nm]).canonical()

    def get_vd(self, Np, Nm=None):
        """Time-domain voltage drop between nodes"""

        return self.get_Vd(Np, Nm).time()
