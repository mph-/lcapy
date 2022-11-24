"""This module provides the NetlistNamespace class.

Copyright 2020--2022 Michael Hayes, UCECE

"""


class NetlistNamespace(object):
    """This class allows elements, nodes, or other namespaces
    in a heirachical namespace to be accessed by name, via __getattr__.

    For example,

    >>> a = Circuit('''
    ... b.L1 1 2'''
    >>> a.b.L1.v

    """

    def __init__(self, namespace, netlist):

        # The namespace name, e.g., a or a.b
        self.namespace = namespace
        self._netlist = netlist
        # The children namespaces
        self.namespaces = {}
        self.allow_anon = False

    def __getitem__(self, name):
        """Return element or node by name."""

        netlist = self._netlist

        # If name is an integer, convert to a string.
        if isinstance(name, int):
            name = '%d' % name

        name = self.namespace + '.' + name

        if name in netlist.nodes:
            return netlist.nodes[name]

        if name in netlist._elements:
            return netlist._elements[name]

        # Try first anonymous name.
        if name + 'anon1' in netlist._elements:
            return netlist._elements[name + 'anon1']

        if name in self.namespaces:
            return self.namespaces[name]

        raise AttributeError('Unknown element or node name %s' % name)

    def __getattr__(self, attr):
        """Return element, node, or another NetlistNamespace object by name.
        This gets called if there is no explicit attribute attr for
        this instance.  This is primarily for accessing elements and
        non-numerical node names.  It also gets called if the called
        attr throws an AttributeError exception.  The annoying thing
        is that hasattr uses getattr and checks for an exception.

        """

        return self.__getitem__(attr)

    @property
    def components(self):
        """Return ordered dictionary of components."""

        return self.elements

    def netlist(self):
        """Return the current netlist for this namespace."""

        nlist = self._netlist

        return '\n'.join([str(cpt) for cpt in nlist._elements.values() if str(cpt).startswith(self.namespace)])

    def __repr__(self):

        return self.netlist()

    @property
    def sch(self):
        """Generate schematic for this namespace."""

        from .schematic import Schematic

        if hasattr(self, '_sch'):
            return self._sch

        sch = Schematic(allow_anon=self.allow_anon)

        netlist = self._netlist.netlist()
        for net in netlist.split('\n'):
            if net.startswith(self.namespace):
                sch.add(net)

        self._sch = sch
        return sch

    def draw(self, filename=None, **kwargs):
        """Draw schematic of subnetlist.

        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.

        Note, if using Jupyter, then need to first issue command %matplotlib inline

        kwargs include:
           label_ids: True to show component ids
           label_values: True to display component values
           draw_nodes:
              True to show all nodes,
              False or 'none' to show no nodes,
              'primary' to show primary nodes,
              'connections' to show nodes that connect more than two components,
              'all' to show all nodes
           label_nodes:
               True to label all nodes,
               False or 'none' to label no nodes,
               'primary' to label primary nodes,
               'alpha' to label nodes starting with a letter,
               'pins' to label nodes that are pins on a chip,
               'all' to label all nodes
           style: 'american', 'british', or 'european'
           scale: schematic scale factor, default 1.0
           node_spacing: spacing between component nodes, default 2.0
           cpt_size: size of a component, default 1.5
           dpi: dots per inch for png files
           help_lines: distance between lines in grid, default 0.0 (disabled)
           debug: True to display debug information
        """

        cct = self
        if kwargs.pop('s_model', False):
            cct = cct._s_model()

        return cct.sch.draw(filename=filename, **kwargs)
