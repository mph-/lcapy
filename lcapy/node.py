"""This module provides support for nodes.

Copyright 2020 Michael Hayes, UCECE

"""

from .immittancemixin import ImmittanceMixin

class Node(ImmittanceMixin):

    def __init__(self, cct, name):

        self.cct = cct
        self.name = name
        self.pos = None
        self.port = False
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        # List of elements connected to this node.
        self._connected = []

    def __repr__(self):
        return "Node('%s')" % self.name
        
    @property
    def V(self):
        """Node voltage with respect to ground."""

        return self.cct.get_Vd(self.name, '0')

    @property
    def v(self):
        """Node time-domain voltage with respect to ground."""

        return self.cct.get_vd(self.name, '0')

    @property
    def admittance(self):
        """Driving-point admittance between node and ground."""

        return self.cct.admittance(self.name, '0')

    @property
    def impedance(self):
        """Driving-point impedance between node and ground."""        

        return self.cct.impedance(self.name, '0')        

    def append(self, cpt):

        if cpt.type in ('P', ):
            self.port = True

        self._connected.append(cpt)

    def oneport(self, node=0):
        """Create oneport object with respect to specified node
        (default ground)."""

        return self.cct.oneport(self.name, node)        
        
    def thevenin(self, node=0):
        """Create Thevenin oneport object with respect to specified
        node (default ground)."""
        
        return self.cct.thevenin(self.name, node)

    def norton(self, node=0):
        """Create Norton oneport object with respect to specified node
        (default ground)."""
        
        return self.cct.norton(self.name, node)    

    @property
    def connected(self):
        """Return list of components connected to the node."""
        
        return self._connected

    def is_connected(self, cpt):
        """Return True if cpt is connected to the node."""

        if isinstance(cpt, str):
            for cpt1 in self.connected:
                if cpt1.name == cpt:
                    return True
            return False
        
        return cpt in self.connected

                
