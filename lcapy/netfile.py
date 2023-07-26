"""This module provides the NetfileMixin class.  This is used for
Netlist and Schematic to parse a netlist file.

Copyright 2020--2023 Michael Hayes, UCECE

"""

from . import grammar
from .parser import Parser
from .state import state
from .componentnamer import ComponentNamer
from warnings import warn
from os.path import dirname, join


class NetfileMixin(object):
    """This parses netlist files for Netlist and Schematic."""

    def _init_parser(self, cpts, allow_anon=False):
        self.parser = Parser(cpts, grammar, allow_anon)
        # Current namespace
        self.namespace = ''
        self.namer = ComponentNamer()
        self.dirname = None
        self.subnetlists = {}

    def _make_anon_cpt_name(self, cpt_type):
        """Make name for anonymous component"""

        return self.namer.name(cpt_type + 'anon', self.elements)

    def _make_anon_node_name(self):
        """Make name for anonymous node"""

        return self.namer.name('nodeanon', self.nodes)

    def _include(self, string):

        parts = string.split(' ')
        if len(parts) < 2 or parts[0] != '.include':
            raise ValueError('Expecting include filename in %s' % string)
        filename = parts[1]
        if len(parts) == 2:
            return self._netfile_add(filename, self.namespace)

        if len(parts) != 4 and parts[2] != 'as':
            raise ValueError(
                'Expecting include filename as name in %s' % string)

        name = parts[3]
        namespace = self.namespace
        self.namespace = name + '.' + namespace
        if name in self.subnetlists:
            warn('Overriding subnetlist %s with %s for %s' %
                 (self.subnetlists[name], parts[1], name))
        self.subnetlists[name] = parts[1]
        ret = self._netfile_add(filename, self.namespace)
        self.namespace = namespace
        return ret

    def _parse(self, string, namespace=''):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.

        This strips leading ...

        `.include filename` will load and parse the specified filename
        `.pdb` enters the debugger
        """

        if string.startswith('...'):
            string = string[3:].strip()

        if string == '':
            pass
        elif string[0:9] == '.include ':
            self._include(string)
            return None
        elif string[0:4] == '.pdb':
            import pdb
            pdb.set_trace()

        cpt = self.parser.parse(string, namespace, self)
        return cpt

    def add(self, string):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Name is the component name, Np is the positive node, and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.

        This returns the added component.
        """

        # Switch context to capture new symbol definitions
        if self.context is not None:
            state.switch_context(self.context)
        cpt = self._add(string)
        self._invalidate()
        if self.context is not None:
            state.restore_context()
        return cpt

    def _add(self, string, namespace=''):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        string = string.strip()
        if '\n' in string:
            lines = string.split('\n')
            for line in lines:
                self._add(line.strip(), namespace)
            return None

        cpt = self._parse(string, namespace)
        if cpt is not None:
            self._cpt_add(cpt)
        return cpt

    def last_added(self):
        """Return name of last added component."""

        return list(self.elements.keys())[-1]

    def netfile_add(self, pathname):
        """Add the nets from file with specified pathname"""

        from os.path import expanduser

        pathname = expanduser(pathname)

        self.dirname = dirname(pathname)
        self._netfile_add(pathname)

    def _netfile_add(self, pathname, namespace=''):
        """Add the nets from file with specified pathname"""

        netfile = None

        try:
            netfile = open(pathname, 'r')
        except:
            pass

        if netfile is None:
            try:
                netfile = open(pathname + '.sch', 'r')
            except:
                pass

        if netfile is None:
            try:
                netfile = open(join(self.dirname, pathname), 'r')
            except:
                pass

        if netfile is None:
            raise FileNotFoundError('Could not open ' + pathname)

        lines = netfile.readlines()
        netfile.close()

        if self.context is not None:
            state.switch_context(self.context)
        for line in lines:
            self._add(line, namespace)
        if self.context is not None:
            state.restore_context()
