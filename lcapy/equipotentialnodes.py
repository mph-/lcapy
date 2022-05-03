"""This module provides the EquipotentialNodes class.  This is used to
find nodes connected by wires.

Copyright 2020--2022 Michael Hayes, UCECE

"""


class EquipotentialNodes(dict):

    def add(self, nodenames):
        for key in nodenames:
            self[key] = [key]

    def find_key(self, n):

        for key, nodes in self.items():
            if n in nodes:
                return key
        return None

    def add_wire(self, n1, n2):

        key1 = self.find_key(n1)
        key2 = self.find_key(n2)

        if key1 is None or key2 is None:
            return

        if key1 != key2:
            # Merge equipotential nodes.  Try to use
            # primary nodes as the key.
            if '_' in key2:
                self[key1].extend(self.pop(key2))
            else:
                self[key2].extend(self.pop(key1))

    def add_wires(self, nodes):

        # Iterate over all pairs.
        for i in range(len(nodes)):
            for m in range(i + 1, len(nodes)):
                self.add_wire(nodes[i], nodes[m])
