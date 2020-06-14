"""
Copyright 2014--2020 Michael Hayes, UCECE
"""

from __future__ import division
from .expr import expr
from .printing import latex, pretty
from .schematic import Schematic
from .circuit import Circuit
from .state import state

class Network(object):
    """This is the base class for network objects."""

    voltage_source = False
    current_source = False

    # True if initial conditions are zero (or unspecified).
    zeroic = True

    # None if component does not have initial conditions.
    # True if initial conditions are specified.
    # False if initial conditions are not specified.
    hasic = None

    netname = ''
    netkeyword = ''

    def _make_id(self, kind):
        """Make identifier"""

        if kind not in self._anon:
            self._anon[kind] = 0
        self._anon[kind] += 1
        return self._anon[kind]
        
    def _tweak_args(self):

        if not hasattr(self, 'args'):
            return ()

        modargs = []
        for arg in self.args:
            arg = expr(arg)
            modargs.append(arg)
        return modargs

    def __repr__(self):

        argsrepr = ', '.join([repr(arg) for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def __str__(self):

        argsrepr = ', '.join([str(arg) for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def _repr_pretty_(self, p, cycle):

        p.text(self.pretty())

    def _repr_latex_(self):

        return '$%s$' % self.latex()

    def pretty(self, **kwargs):

        argsrepr = ', '.join([pretty(arg, **kwargs) for arg in self._tweak_args()])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def latex(self, **kwargs):

        argsrepr = ', '.join([latex(arg, **kwargs) for arg in self._tweak_args()])
        return '\\mathrm{%s}(%s)' % (self.__class__.__name__, argsrepr)

    @property
    def analysis(self):
        return self.cct.analysis

    def describe(self):
        """Print a message describing how network is solved."""
        return self.cct.describe()
    
    def simplify(self):

        return self

    def _add_elements(self):

        netlist = self.netlist()
        for net in netlist.split('\n'):
            self._add(net)

        # Hack, create ground reference.
        self._add('W %d 0' % (self.node - 1))

    @property 
    def node(self):

        if not hasattr(self, 'node_counter'):
            self.node_counter = 0
        ret = self.node_counter
        self.node_counter += 1
        return ret

    def netargs(self):

        def quote(arg):

            # TODO: make more robust to catch expressions.
            if ('(' in arg) or (')' in arg) or (' ' in arg) or (',' in arg) or ('*' in arg) or ('/' in arg):
                return '{%s}' % arg
            return arg

        return ' '.join([quote(str(arg)) for arg in self.args])

    def net_make(self, net, n1=None, n2=None):

        if n1 == None:
            n1 = net.node
        if n2 == None:
            n2 = net.node

        netname = self.__class__.__name__ if self.netname == '' else self.netname

        netid = net._make_id(netname)
        if self.netkeyword != '':
            return '%s%s %s %s %s %s; right' % (netname, netid,
                                                n1, n2, 
                                                self.netkeyword, self.netargs())
        else:
            return '%s%s %s %s %s; right' % (netname, netid,
                                             n1, n2, self.netargs())

    def netlist(self):

        # Enumerate from node 0
        self.node_counter = 0
        self._anon = {}
        n1 = self.node
        n2 = self.node        
        return self.net_make(self, n2, n1)

    def pdb(self):
        """Enter the python debugger."""
        
        import pdb; pdb.set_trace()
        return self
    
    @property
    def sch(self):
        """Convert a Network object into a Schematic object."""

        if hasattr(self, '_sch'):
            return self._sch

        netlist = self.netlist()
        sch = Schematic()
        for net in netlist.split('\n'):
            sch.add(net)
        self._sch = sch
        return sch

    def draw(self, filename=None, **kwargs):
        """Draw schematic of network.

        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.

        Note, if using Jupyter, then need to first issue command %matplotlib inline

        kwargs include:
           label_ids: True to show component ids
           label_values: True to display component values
           draw_nodes: True to show all nodes, False to show no nodes, 
             'primary' to show primary nodes,
             'connections' to show nodes that connect more than two components,
             'all' to show all nodes
           label_nodes: True to label all nodes, False to label no nodes, 
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

        if 'label_ids' not in kwargs:
            kwargs['label_ids'] = False
        if 'label_values' not in kwargs:
            kwargs['label_values'] = True
        if 'label_nodes' not in kwargs:
            kwargs['label_nodes'] = False
        if 'draw_nodes' not in kwargs:
            kwargs['draw_nodes'] = 'connections'
        
        self.sch.draw(filename=filename, **kwargs)
        
    @property
    def cct(self):
        """Convert a Network object into a Circuit object."""

        if hasattr(self, '_cct'):
            return self._cct

        netlist = self.netlist()
        cct = Circuit()
        for net in netlist.split('\n'):
            cct.add(net)
        self._cct = cct
        return cct

    def circuit(self):
        """Convert a Network object into a Circuit object."""

        return self.cct

    @property    
    def initial_value_problem(self):
        return self.cct.initial_value_problem

    @property    
    def is_ivp(self):
        return self.cct.is_ivp

    @property
    def is_dc(self):
        return self.cct.is_dc

    @property    
    def is_ac(self):
        return self.cct.is_ac

    @property    
    def is_causal(self):
        return self.cct.is_causal

    @property
    def has_dc(self):
        return self.cct.has_dc

    @property    
    def has_ac(self):
        return self.cct.has_ac

    @property    
    def has_transient(self):
        return self.cct.has_transient

    @property    
    def kinds(self):
        """Return list of transform domain kinds."""        
        return self.cct.kinds

    @property
    def symbols(self):
        """Return dictionary of symbols defined in the network."""
        
        return self.cct.symbols

    @property
    def all_symbols(self):
        """Return dictionary of symbols defined in the network and global
        symbols."""

        symbols = self.symbols
        symbols.update(state.global_context.symbols)
        return symbols
