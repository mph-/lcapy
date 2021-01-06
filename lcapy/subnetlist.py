"""This module provides the SubNetlist class.

Copyright 2020 Michael Hayes, UCECE

"""

from .mna import MNAMixin
from .netlistmixin import NetlistMixin
from .state import state
from .symbols import omega



class SubNetlist(NetlistMixin, MNAMixin):
    """This is a representation of a netlist for a particular
    transformation domain, such as ac, dc, transient, or noise."""

    def __new__(cls, netlist, kind):

        obj = netlist.select(kind=kind)
        # Need own context to avoid conflicts with Vn1 and Vn1(s), etc.
        obj.context = state.new_context()
        obj.kind = kind
        obj.__class__ = cls
        obj._analysis = obj.analyse()
        return obj

    def __init__(self, netlist, kind):
        """ kind can be 't', 'dc', 's', 'time', 'ivp', 'n*' or omega, 
        where 'n*' is a noise identifer and omega is an angular frequency."""
        
        if not isinstance(kind, str):
            return
        if kind[0] == 'n':
            return
        kinds = ('t', 'dc', 's', 'time', 'ivp', 'laplace')
        if kind not in kinds:
            raise ValueError('Expected one of %s for kind, got %s' %
                             (', '.join(kinds), kind))

    def get_I(self, name):
        """Current through component"""

        self._solve()
        return self._Idict[name].canonical()

    def get_i(self, name):
        """Time-domain current through component"""

        return self.get_I(name).time()

    def get_Vd(self, Np, Nm=None):
        """Voltage drop between nodes"""

        self._solve()
        return (self._Vdict[Np] - self._Vdict[Nm]).canonical()

    def get_vd(self, Np, Nm=None):
        """Time-domain voltage drop between nodes"""

        return self.get_Vd(Np, Nm).time()
