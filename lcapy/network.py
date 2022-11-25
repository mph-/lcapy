"""
Copyright 2014--2022 Michael Hayes, UCECE
"""

from __future__ import division
from .expr import expr
from .sexpr import s
from .printing import latex, pretty
from .circuit import Circuit
from .state import state
from warnings import warn


class Network(object):
    """This is the base class for network objects."""

    is_voltage_source = False
    is_current_source = False
    is_inductor = False
    is_capacitor = False
    is_resistor = False
    is_conductor = False
    is_parallel = False
    is_series = False
    is_noiseless = True

    # True if initial conditions are zero (or unspecified).
    zeroic = True

    # None if component does not have initial conditions.
    # True if initial conditions are specified.
    # False if initial conditions are not specified.
    has_ic = None

    cpt_type = ''
    netkeyword = ''
    kwargs = {}

    def _opts_str(self, arg=None, **kwargs):

        parts = []
        if arg is not None and arg != '':
            parts.append(arg)

        for key, val in kwargs.items():
            parts.append('%s=%s' % (key, val))

        for key, val in self.kwargs.items():
            parts.append('%s=%s' % (key, val))
        return ', '.join(parts)

    def __repr__(self):

        parts = []
        for key, val in self.kwargs.items():
            parts.append("%s='%s'" % (key, val))

        argsrepr = ', '.join([repr(arg) for arg in self.args] + parts)
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def __str__(self):

        argsrepr = ', '.join([repr(arg) for arg in self.args])
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def _repr_pretty_(self, p, cycle):

        p.text(self.pretty())

    def _repr_latex_(self):

        s = self.latex(mode='plain')
        return "$$%s$$" % s

    def _repr_png_(self):

        return None

    def _repr_svg_(self):

        return None

    def pretty(self, **kwargs):

        def pretty1(arg, **kwargs):
            try:
                return arg.pretty(**kwargs)
            except:
                return pretty(arg, **kwargs)

        parts = []
        for key, val in self.kwargs.items():
            parts.append("%s='%s'" % (key, val))

        argsrepr = ', '.join([pretty1(arg, **kwargs)
                             for arg in self.args] + parts)
        return '%s(%s)' % (self.__class__.__name__, argsrepr)

    def latex(self, **kwargs):

        def latex1(arg, **kwargs):
            try:
                return arg.latex(**kwargs)
            except:
                return latex(arg, **kwargs)

        argsrepr = ', '.join([latex1(arg, **kwargs) for arg in self.args])
        return '\\mathrm{%s}(%s)' % (self.__class__.__name__, argsrepr)

    def pprint(self):

        print(self.pretty())

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
        self._add('W %d 0' % (self._node - 1))

    def _net_make(self, netlistmaker, n1=None, n2=None, dir='right'):

        net = self
        if n2 == None:
            n2 = net._node
        if n1 == None:
            n1 = net._node

        cpt_type = net.__class__.__name__ if net.cpt_type == '' else net.cpt_type

        opts_str = self._opts_str(dir)
        name = netlistmaker._make_name(cpt_type, self.args)
        if net.netkeyword != '':
            return '%s %s %s %s %s; %s' % (name, n1, n2, net.netkeyword,
                                           netlistmaker._netargs(net), opts_str)
        else:
            return '%s %s %s %s; %s' % (name, n1, n2,
                                        netlistmaker._netargs(net), opts_str)

    @property
    def _depths(self):
        return [net._depth for net in self.args]

    @property
    def _depth(self):
        from .oneport import Ser, Par

        if not isinstance(self, (Ser, Par)):
            return 0

        depths = self._depths
        return 1 + max(depths)

    def netlist(self, layout='horizontal', evalf=None):
        """Create a netlist.

        `layout` can be 'horizontal', 'vertical', or 'ladder'.

        `evalf` can be False or an integer specifying the number of
        decimal places used to evaluate floats.
        """

        if layout == 'ladder':
            from .laddermaker import LadderMaker
            return LadderMaker(self, layout=layout, evalf=evalf)()

        from .netlistmaker import NetlistMaker
        return NetlistMaker(self, layout=layout, evalf=evalf)()

    def pdb(self):
        """Enter the python debugger."""

        import pdb
        pdb.set_trace()
        return self

    def sch(self, layout='horizontal', evalf=False):
        """Convert a Network object into a Schematic object."""

        from .schematic import Schematic

        netlist = self.netlist(layout=layout, evalf=evalf)
        sch = Schematic()
        for net in netlist.split('\n'):
            sch.add(net)

        return sch

    def draw(self, filename=None, layout='horizontal', form=None, evalf=False, **kwargs):
        """Draw schematic of network.

        filename specifies the name of the file to produce.  If None,
        the schematic is displayed on the screen.

        Note, if using Jupyter, then need to first issue command %matplotlib inline

        `layout` is either 'horizontal', 'vertical', or 'ladder'.

        `evalf` can be False or an integer specifying the number of
        decimal places used to evaluate floats.

        `kwargs` include:
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

        if form is not None:
            warn('Form is deprecated; use layout instead')
            layout = form

        self.sch(layout=layout, evalf=evalf).draw(filename=filename, **kwargs)

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
    def is_IVP(self):
        return self.cct.is_IVP

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
    def params(self):
        """Return list of symbols used as arguments in the network."""

        # This would be faster to determine directly from AST.
        return self.cct.params

    @property
    def symbols(self):
        """Return dictionary of symbols defined in the network."""

        return self.cct.symbols

    def transform(self, form='cauerI'):
        """Transform the network impedance into an alternative form.  The
        transformation is performed using network synthesis of the
        network's impedance (note, this ignores the sources).  `form`
        includes: cauerI, cauerII, fosterI, fosterII.

        Note some methods generate networks with negative value
        components.

        """

        return self.Z(s).network(form)

    def subs(self, subs_dict):
        """Substitute values using dictionary of substitutions.

        For example, b = a.subs({'R1': 1e3, 'R2': 9e3})"""

        if not self.is_parallel and not self.is_series:
            return self.__class__(*[str(expr(arg).subs(subs_dict)) for arg in self.args])

        return self.__class__(*[arg.subs(subs_dict) for arg in self.args])

    def noisy(self, T='T'):
        """Create noisy network model by replacing resistances with a series
        combination of a resistance and a noise voltage source."""

        if self.is_resistor and not self.is_noiseless:
            from .oneport import Vnoise, NR

            Vn = 'sqrt(4 * k_B * %s * %s)' % (T, self.args[0])
            return NR(*self.args) + Vnoise(Vn)

        if self.is_conductor and not self.is_noiseless:
            from .oneport import Vnoise, NG

            Vn = 'sqrt(4 * k_B * %s / (%s))' % (T, self.args[0])
            return NG(*self.args) + Vnoise(Vn)

        if not self.is_parallel and not self.is_series:
            return self.__class__(*self.args)

        return self.__class__(*[arg.noisy(T=T) for arg in self.args])
