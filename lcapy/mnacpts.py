
"""
This module defines the components for modified nodal analysis.  The components
are defined at the bottom of this file.

Copyright 2015--2021 Michael Hayes, UCECE

"""

from __future__ import print_function
from .expr import expr
from .cexpr import ConstantDomainExpression
from .omegaexpr import AngularFourierDomainExpression
from .symbols import j, omega, jomega, s, t
from .functions import sqrt
from .sym import capitalize_name, omegasym
from .grammar import delimiters
from .immittancemixin import ImmittanceMixin
from .superpositioncurrent import SuperpositionCurrent
from .voltage import voltage
from .current import current
from .opts import Opts
import lcapy
import inspect
import sys
import sympy as sym

module = sys.modules[__name__]

cptaliases = {'E': 'VCVS', 'F': 'CCCS',
              'G': 'VCCS',  'H': 'CCVS',
              'r': 'Damper', 'm': 'Mass',
              'k': 'Spring'}


class Cpt(ImmittanceMixin):

    dependent_source = False
    independent_source = False    
    reactive = False
    need_branch_current = False
    need_extra_branch_current = False    
    need_control_current = False
    directive = False
    flip_branch_current = False
    ignore = False
    add_series = False
    add_parallel = False    
    equipotential_nodes = ()
    

    def __init__(self, cct, namespace, defname, name, cpt_type, cpt_id, string,
                 opts_string, nodes, keyword, *args):

        self.cct = cct
        self.type = cpt_type
        self.id = cpt_id
        self.defname = defname
        self.name = name
        self.relname = name
        self.namespace = ''
        self.nodenames = nodes
        self.relnodes = nodes
        
        parts = name.split('.')
        if len(parts) > 1:
            self.namespace = '.'.join(parts[0:-1]) + '.'
            self.relname = parts[-1]
            self.relnodes = []
            for node in nodes:
                if node.startswith(self.namespace):
                    node = node[len(self.namespace):]
                self.relnodes.append(node)

        self._string = string
        #self.net = string.split(';')[0]
        # This is the initial opts_string from which the opts attribute
        # is derived.
        self.opts_string = opts_string
        self.args = args
        self.explicit_args = args        
        self.classname = self.__class__.__name__
        self.keyword = keyword
        self.opts = Opts(self.opts_string)

        # No defined cpt
        if self.type in ('XX', 'Cable'):
            self._cpt = lcapy.oneport.Dummy()
            return

        if ((args == () and not self.type in ('W', 'O', 'P'))
            or (self.type in ('F', 'H', 'CCCS', 'CCVS') and len(args) == 1)
            or (self.type == 'K' and len(args) == 2)):
            # Default value is the component name
            value = self.type + self.id

            if False and self.type in ('V', 'I') and keyword[1] == '':
                # This is too subtle and creates havoc with
                # the expected behaviour of subs.
                value = value[0].lower() + value[1:] + '(t)'

            args += (value, )
            self.args = args

        classname = self.classname
        # Handle aliases.
        try:
            classname = cptaliases[classname]
        except:
            pass
                
        try:
            newclass = getattr(lcapy.oneport, classname)
        except:
            try:
                newclass = getattr(lcapy.twoport, classname)
            except:
                self._cpt = lcapy.oneport.Dummy()
                return
                
        self._cpt = newclass(*args)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        # TODO, use self._netmake() but fix for XX
        return self._string        

    @property
    def cpt(self):
        return self._cpt
    
    def _stamp(self, cct):
        raise NotImplementedError('stamp method not implemented for %s' % self)

    def _copy(self):
        """Make copy of net."""
        
        return str(self)
    
    def _arg_format(self, value):
        """Place value string inside curly braces if it contains a delimiter."""

        string = str(value)

        if string.startswith('{'):
            return string
        
        for delimiter in delimiters:
            if delimiter in string:
                return '{' + string + '}'
        return string

    def _kill(self):
        """Kill sources."""

        return self._copy()

    def _subs(self, subs_dict):
        """Substitute values using dictionary of substitutions.
        If a scalar is passed, this is substituted for the component value.

        For example, given a component, cpt, defined as 'R1 1 2' then
        cpt.subs(5) and cpt.subs({'R1': 5}) are equivalent.  In both
        cases, the result is 'R1 1 2 5'."""

        if not isinstance(subs_dict, dict):
            subs_dict = {self.args[0]: subs_dict}

        return self._netsubs(subs_dict=subs_dict)

    def _initialize(self, ic):
        """Change initial condition to ic."""

        return self._copy()    

    def _netsubs(self, node_map=None, zero=False, subs_dict=None):
        """Create a new net description.  If node_map is not None,
        rename the nodes.  If zero is True, set args to zero."""

        string = self.defname
        field = 0
        
        for node in self.relnodes:
            if node_map is not None:
                node = node_map[node]
            string += ' ' + node
            field += 1
            if field == self.keyword[0]:
                string += ' ' + self.keyword[1]
                field += 1

        if subs_dict is None:
            args = self.explicit_args
        else:
            args = self.args
                
        for arg in self.args:
            if zero:
                arg = 0

            # Perform substitutions
            if subs_dict is not None:
                arg = str(expr(arg).subs(subs_dict))
                
            string += ' ' + self._arg_format(arg)
            field += 1
            if field == self.keyword[0]:
                string += self.keyword[1]
        opts_str = str(self.opts).strip()
        if opts_str != '':
            string += '; ' + opts_str
        return string
        
    def _rename_nodes(self, node_map):
        """Rename the nodes using dictionary node_map."""

        return self._netsubs(node_map)

    def _netmake1(self, name, nodes=None, args=None, opts=None):

        if nodes is None:
            nodes = self.relnodes

        if args is None:
            args = self.args

        if not isinstance(args, tuple):
            args = (args, )

        if opts is None:
            opts = self.opts
            
        fmtargs = []
        for arg in args:
            fmtargs.append(self._arg_format(arg))

        if len(fmtargs) == 1 and fmtargs[0] == name:
            fmtargs = []

        parts = [name]
        parts.extend(nodes)
        parts.extend(fmtargs)

        # Insert keyword...
        
        
        net = ' '.join(parts)
        opts_str = str(opts).strip()
        if opts_str != '':
            net += '; ' + opts_str
        return net
    
    def _netmake(self, nodes=None, args=None, opts=None):
        """This keeps the same cpt name"""                
        return self._netmake1(self.namespace + self.relname, nodes, args, opts)

    def _netmake_W(self, nodes=None, opts=None):
        """This is used for changing cpt name from L1 to W"""        
        return self._netmake1(self.namespace + 'W', nodes, args=(),
                              opts=opts)

    def _netmake_O(self, nodes=None, opts=None):
        """This is used for changing cpt name from C1 to O"""        
        return self._netmake1(self.namespace + 'O', nodes, args=(),
                              opts=opts)    
        
    def _netmake_variant(self, newtype, nodes=None, args=None, opts=None,
                         suffix=''):
        """This is used for changing cpt name from C1 to ZC1"""

        name = self.namespace + newtype + self.relname + suffix
        return self._netmake1(name, nodes, args, opts)
    
    def _select(self, kind=None):
        """Select domain kind for component."""
        raise ValueError('Component not a source: %s' % self)    

    def _new_value(self, value, ic=None):

        args = (value, )
        if ic is not None:
            args = (value, ic)
        
        return self._netmake(args=args)
    
    def _zero(self):
        """Zero value of the voltage source.  This kills it but keeps it as a
        voltage source in the netlist.  This is required for dummy
        voltage sources that are required to specify the controlling
        current for CCVS and CCCS components."""        

        raise ValueError('Component not a source: %s' % self)        

    def _r_model(self):
        """Return resistive model of component."""
        return self._copy()

    def _s_model(self, var):
        """Return s-domain model of component."""
        return self._copy()    

    def _ss_model(self):
        """Return state-space model of component."""
        return self._copy()            

    def _pre_initial_model(self):
        """Return pre-initial model of component."""
        return self._copy()        

    @property
    def is_source(self):
        """Return True if component is a source (dependent or independent)"""
        return self.dependent_source or self.independent_source

    @property
    def is_dependent_source(self):
        """Return True if component is a dependent source"""
        return self.dependent_source

    @property
    def is_independent_source(self):
        """Return True if component is an independent source"""
        return self.independent_source

    @property
    def _source_IV(self):

        if self.cpt.is_voltage_source:
            return self.cpt.Voc
        elif self.cpt.is_current_source:
            return self.cpt.Isc
        else:
            raise ValueError('%s is not a source' % self)        

    @property
    def is_causal(self):
        """Return True if causal component or if source produces
        a causal signal."""

        if self.is_source:
            return self._source_IV.is_causal
        return self.cpt.is_causal

    @property
    def is_dc(self):
        """Return True if source is dc."""

        return self._source_IV.is_dc

    @property
    def is_ac(self):
        """Return True if source is ac."""

        return self._source_IV.is_ac

    @property
    def has_ac(self):
        """Return True if source has ac component."""

        return self._source_IV.has_ac

    @property
    def has_dc(self):
        """Return True if source has dc component."""

        return self._source_IV.has_dc    

    @property
    def has_noisy(self):
        """Return True if source has noisy component."""

        return self._source_IV.has_noisy
    
    @property
    def has_s_transient(self):
        """Return True if source has transient component defined in s-domain."""

        return self._source_IV.has_s_transient

    @property
    def has_t_transient(self):
        """Return True if source has transient component defined in time domain."""

        return self._source_IV.has_t_transient        

    @property
    def has_transient(self):
        """Return True if source has a transient component."""

        return self._source_IV.has_transient        
        
    @property
    def is_noisy(self):
        """Return True if source is noisy."""

        if self.cpt.is_voltage_source:
            return self.cpt.is_noisy
        elif self.cpt.is_current_source:
            return self.cpt.is_noisy
        else:
            raise ValueError('%s is not a source' % self)        

    @property
    def is_noiseless(self):
        """Return True if component is noiseless."""

        return self.cpt.is_noiseless
        
    @property
    def is_inductor(self):
        """Return True if component is an inductor."""
        return self.cpt.is_inductor

    @property
    def is_capacitor(self):
        """Return True if component is a capacitor."""
        return self.cpt.is_capacitor

    @property
    def is_reactance(self):
        """Return True if component is a capacitor or inductor."""
        return self.is_capacitor or self.is_inductor

    @property
    def is_resistor(self):
        """Return True if component is a resistor."""
        return self.cpt.is_resistor

    @property
    def is_conductor(self):
        """Return True if component is a conductor."""
        return self.cpt.is_conductor            

    @property
    def is_voltage_source(self):
        """Return True if component is a voltage source (dependent or
        independent)"""
        return self.cpt.is_voltage_source    

    @property
    def is_current_source(self):
        """Return True if component is a current source (dependent or
        independent)"""
        return self.cpt.is_current_source    
        
    @property
    def zeroic(self):
        """Return True if initial conditions are zero (or unspecified)."""

        return self.cpt.zeroic

    @property
    def has_ic(self):
        """Return True if initial conditions are specified."""

        return self.cpt.has_ic

    @property
    def I(self):
        """Current through component.  The current is defined to be into the
        positive node for passive devices and out of the positive node
        for sources."""

        return self.cct.get_I(self.name)

    @property
    def i(self):
        """Time-domain current through component.  The current is
        defined to be into the positive node for passive devices and
        out of the positive node for sources."""

        return self.cct.get_i(self.name)

    @property
    def V(self):
        """Voltage drop across component."""

        return self.cct.get_Vd(self.nodenames[0], self.nodenames[1])

    @property
    def v(self):
        """Time-domain voltage drop across component."""

        return self.cct.get_vd(self.nodenames[0], self.nodenames[1])

    @property
    def Isc(self):
        """Short-circuit current for component in isolation, i.e, current in
        wire connected across component."""

        return self.cpt.Isc.select(self.cct.kind)

    @property
    def isc(self):
        """Short-circuit time-domain current for component in isolation, i.e,
        current in wire connected across component."""

        return self.Isc.time()

    @property
    def Voc(self):
        """Open-circuit voltage for component in isolation."""

        return self.cpt.Voc.select(self.cct.kind)

    @property
    def voc(self):
        """Open-circuit time-domain voltage for component in isolation."""

        return self.Voc.time()    
    
    @property
    def V0(self):
        """Initial voltage (for capacitors only)."""

        return voltage(0)

    @property
    def I0(self):
        """Initial current (for inductors only)."""

        return current(0)
    
    @property
    def admittance(self):
        """Self admittance of component. 

        The admittance is expressed in jomega form for AC circuits
        and in s-domain for for transient circuits.

        For the driving-point admittance measured across the component
        use .dpY or .oneport().Y"""

        return self.cpt.admittance._select(self.cct.kind)        

    @property
    def impedance(self):
        """Self impedance of component.  

        The impedance is expressed in jomega form for AC circuits
        and in s-domain for for transient circuits.

        For the driving-point impedance measured across the component
        use .dpZ or .oneport().Z"""        

        return self.cpt.impedance._select(self.cct.kind)

    @property
    def dpIsc(self):
        """Driving-point short-circuit current, i.e, current in wire
        connected across component connected in-circuit.

        """
        return self.oneport().Isc

    @property
    def dpisc(self):
        """Driving-point short-circuit time-domain current i.e, current in
        wire connected across component in-circuit."""

        return self.dpIsc.time()

    @property
    def dpVoc(self):
        """Driving-point open-circuit voltage across component in-circuit."""

        return self.oneport().V

    @property
    def dpvoc(self):
        """Driving-point open-circuit time-domain voltage across component in
        circuit."""

        return self.dpVoc.time()    

    @property
    def dpYs(self):
        """Driving-point generalized admittance (s-domain) measured across
        component in-circuit.  For the admittance of the component in isolation
        use .Ys"""

        return self.oneport().Ys

    @property
    def dpZs(self):
        """Driving-point generalized impedance (s-domain) measured across
        component in-circuit.  For the impedance of the component in isolation
        use .Zs"""        

        return self.oneport().Zs
    
    @property
    def dpY(self):
        """Driving-point admittance measured across component in-circuit.  For
        the admittance of the component in isolation use .Y"""        

        return self.oneport().Y

    @property
    def dpZ(self):
        """Driving-point impedance measured across component in-circuit.  For
        the impedance of the component in isolation use .Z"""

        return self.oneport().Z

    @property
    def node_indexes(self):

        return [self.cct._node_index(n) for n in self.nodenames]

    @property
    def branch_index(self):

        return self.cct._branch_index(self.name)

    def dummy_node(self):

        return '_' + self.cct._make_anon('node')

    def oneport(self):
        """Create oneport object."""
        
        return self.cct.oneport(self.nodenames[1], self.nodenames[0])

    def thevenin(self):
        """Create Thevenin oneport object."""
        
        return self.cct.thevenin(self.nodenames[1], self.nodenames[0])

    def norton(self):
        """Create Norton oneport object."""
        
        return self.cct.norton(self.nodenames[1], self.nodenames[0])

    def transfer(self, cpt):
        """Create transfer function for the s-domain voltage across the
        specified cpt divided by the s-domain voltage across self."""

        if isinstance(cpt, str):
            cpt = self.cct._elements[cpt]        
        
        return self.cct.transfer(self.nodenames[1], self.nodenames[0],
                                 cpt.nodenames[1], cpt.nodenames[0])

    @property
    def nodes(self):
        """Return list of nodes for this component.   See also
        nodenames."""

        return [self.cct.nodes[nodename] for nodename in self.nodenames]
    
    def connected(self):
        """Return list of components connected to this component."""

        cpts = set()
        for node in self.nodes:
            cpts = cpts.union(set(node.connected))
        
        return list(cpts)

    def is_connected(self, cpt):
        """Return True if cpt is connected to this component."""

        if isinstance(cpt, str):
            for cpt1 in self.connected:
                if cpt1.name == cpt:
                    return True
            return False
        
        return cpt in self.connected

    
class Invalid(Cpt):
    
    @property
    def cpt(self):
         raise NotImplementedError('Invalid component for circuit analysis: %s' % self)       


class NonLinear(Invalid):

    def _stamp(self, cct):
        raise NotImplementedError('Cannot analyse non-linear component: %s' % self)


class TimeVarying(Invalid):

    def _stamp(self, cct):
        raise NotImplementedError('Cannot analyse time-varying component: %s' % self)


class Logic(Invalid):

    def _stamp(self, cct):
        raise NotImplementedError('Cannot analyse logic component: %s' % self)


class Misc(Invalid):

    def _stamp(self, cct):
        raise NotImplementedError('Cannot analyse misc component: %s' % self)


class Ignored(Cpt):

    ignore = True

    def _stamp(self, cct):
        # Could print warning
        pass
    

class Dummy(Cpt):

    causal = True
    dc = False
    ac = False
    zeroic = True
    has_ic = None
    noisy = False


class XX(Dummy):
    directive = True
    ignore = True
    
    def _stamp(self, cct):
        pass

    def _subs(self, subs_dict):
        return self._copy()

    def _rename_nodes(self, node_map):
        """Rename the nodes using dictionary node_map."""

        return self._copy()

    def __str__(self):        
        return self._string
    
    
class IndependentSource(Cpt):

    independent_source = True
    
    def _zero(self):
        """Zero value of the source.  For a voltage source this makes it a
        short-circuit; for a current source this makes it
        open-circuit.  This effectively kills the source but keeps it
        as a source in the netlist.  This is required for dummy
        voltage sources that are required to specify the controlling
        current for CCVS and CCCS components.

        """
        return self._netsubs(zero=True)
    

class DependentSource(Dummy):

    dependent_source = True        

    def _zero(self):
        return self._copy()

    
class RLC(Cpt):

    def _s_model(self, var):

        if self.Voc == 0:
            return self._netmake_variant('Z', args=self.Z(var))

        dummy_node = self.dummy_node()

        opts = self.opts.copy()

        # Strip voltage labels and save for open-circuit cpt
        # in parallel with Z and V.
        voltage_opts = opts.strip_voltage_labels()

        znet = self._netmake_variant('Z', nodes=(self.relnodes[0], dummy_node),
                                     args=self.Z(var), opts=opts)

        # Strip voltage and current labels from voltage source.
        opts.strip_all_labels()

        vnet = self._netmake_variant('V', nodes=(dummy_node, self.relnodes[1]),
                                     args=self.Voc.laplace()(var), opts=opts)
        if voltage_opts == {}:
            return znet + '\n' + vnet

        # Create open circuit in parallel to the Z and V
        # that has the voltage labels.
        opts = self.opts.copy()
        opts.strip_current_labels()
        # Need to convert voltage labels to s-domain.
        # v(t) -> V(s)
        # v_C -> V_C
        # v_L(t) -> V_L(s)
        for opt, val in voltage_opts.items():
            opts[opt] = capitalize_name(val)

        onet = self._netmake_O(opts=opts)
        return znet + '\n' + vnet + '\n' + onet


class RC(RLC):

    def _noisy(self, T='T'):

        dummy_node = self.dummy_node()

        opts = self.opts.copy()

        # Should check for namespace conflict if user has defined
        # a noiseless resistor.
        rnet = self._netmake_variant('N', nodes=(self.relnodes[0], dummy_node),
                                     args=str(self.R), opts=opts)

        # Use k_B for Boltzmann's constant to avoid clash with k symbol
        # for discrete frequency
        Vn = 'sqrt(4 * k_B * %s * %s)' % (T, self.args[0])
        vnet = self._netmake_variant('Vn', nodes=(dummy_node, self.relnodes[1]),
                                     args=('noise', Vn), opts=opts)
        return rnet + '\n' + vnet
    
    def _stamp(self, cct):

        # L's can also be added with this stamp but if have coupling
        # it is easier to generate a stamp that requires the branch current
        # through the L.
        n1, n2 = self.node_indexes

        if self.type == 'C' and cct.kind == 'dc':
            Y = 0
        else:
            Y = self.Y.expr

        if n1 >= 0 and n2 >= 0:
            cct._G[n1, n2] -= Y
            cct._G[n2, n1] -= Y
        if n1 >= 0:
            cct._G[n1, n1] += Y
        if n2 >= 0:
            cct._G[n2, n2] += Y

        if cct.kind == 'ivp' and self.cpt.has_ic:
            I = self.Isc.expr
            if n1 >= 0:
                cct._Is[n1] += I 
            if n2 >= 0:
                cct._Is[n2] -= I               


class C(RC):

    reactive = True
    add_parallel = True

    @property
    def C(self):
        return self.cpt.C
    
    def _kill(self):
        """Kill implicit sources due to initial conditions."""
        return self.netmake(args=self.args[0])

    def _initialize(self, ic):
        """Change initial condition to ic."""
        return self._netmake(args=(self.args[0], ic))

    def _pre_initial_model(self):

        return self._netmake_variant('V', args=self.cpt.v0)

    def _r_model(self):

        dummy_node = self.dummy_node()
        opts = self.opts.copy()        

        # Use Thevenin model.  This will require the current through
        # the voltage source to be explicitly computed.
        
        Req = 'R%seq' % self.name
        Veq = 'V%seq' % self.name  

        opts.strip_voltage_labels()
        rnet = self._netmake_variant('R', suffix='eq',
                                     nodes=(self.relnodes[0], dummy_node),
                                     args=Req, opts=opts)

        opts.strip_current_labels()        
        vnet = self._netmake_variant('V', suffix='eq',
                                     nodes=(dummy_node, self.relnodes[1]),
                                     args=('dc', Veq), opts=opts)

        # TODO: the voltage labels should be added across an
        # open-circuit object.
        
        return rnet + '\n' + vnet

    def _ss_model(self):
        # Perhaps mangle name to ensure it does not conflict
        # with another voltage source?
        return self._netmake_variant('V_', args='v_%s(t)' % self.relname)        
        
    @property
    def V0(self):
        """Initial voltage (for capacitors only)."""

        if self.cct.kind == 'ivp' and self.cpt.has_ic:
            return voltage(self.cpt.v0 / s)
        return voltage(0)


class CPE(RC):

    # If n == 0, then not reactive
    reactive = True

    pass


class VCVS(DependentSource):
    """VCVS"""

    need_branch_current = True

    def _stamp(self, cct):
        n1, n2, n3, n4 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] += 1
            cct._C[m, n1] += 1
        if n2 >= 0:
            cct._B[n2, m] -= 1
            cct._C[m, n2] -= 1

        Ad = ConstantDomainExpression(self.args[0]).expr
        if len(self.args) > 1:
            Ac = ConstantDomainExpression(self.args[1]).expr
        else:
            Ac = 0

        Ap = (Ac / 2 + Ad)
        Am = (Ac / 2 - Ad)        
            
        if n3 >= 0:
            cct._C[m, n3] -= Ap
        if n4 >= 0:
            cct._C[m, n4] -= Am

    def _kill(self):
        newopts = self.opts.copy()
        newopts.strip_current_labels()
        newopts.strip_labels()

        return self._netmake_W(opts=newopts)

    
class CCCS(DependentSource):
    """CCCS"""

    need_control_current = True
    
    def _stamp(self, cct):
        n1, n2 = self.node_indexes
        m = cct._branch_index(self.args[0])
        F = ConstantDomainExpression(self.args[1]).expr
            
        if n1 >= 0:
            cct._B[n1, m] -= F
        if n2 >= 0:
            cct._B[n2, m] += F

    def _kill(self):
        newopts = self.opts.copy()
        newopts.strip_voltage_labels()
        newopts.strip_labels()

        return self._netmake_O(opts=newopts)


class FB(Misc):
    """Ferrite bead"""
    pass


class VCCS(DependentSource):
    """VCCS"""

    def _stamp(self, cct):
        n1, n2, n3, n4 = self.node_indexes
        G = ConstantDomainExpression(self.args[0]).expr

        if n1 >= 0 and n3 >= 0:
            cct._G[n1, n3] -= G
        if n1 >= 0 and n4 >= 0:
            cct._G[n1, n4] += G
        if n2 >= 0 and n3 >= 0:
            cct._G[n2, n3] += G
        if n2 >= 0 and n4 >= 0:
            cct._G[n2, n4] -= G

    def _kill(self):
        newopts = self.opts.copy()
        newopts.strip_voltage_labels()
        newopts.strip_labels()

        return self._netmake_O(opts=newopts)        

    
class GY(Dummy):
    """Gyrator"""    

    need_branch_current = True
    need_extra_branch_current = True

    def _stamp(self, cct):
        
        n1, n2, n3, n4 = self.node_indexes
        m1 = self.cct._branch_index(self.defname + 'X')
        m2 = self.branch_index

        # m1 is the input branch
        # m2 is the output branch
        # GY.I gives the current through the output branch

        # Could generalise to have different input and output
        # impedances, Z1 and Z2, but if Z1 != Z2 then the device is
        # not passive.

        # V2 = -I1 Z2     V1 = I2 Z1
        # where V2 = V[n1] - V[n2] and V1 = V[n3] - V[n4]
        
        Z1 = ConstantDomainExpression(self.args[0]).expr                    
        Z2 = Z1
        
        if n1 >= 0:
            cct._B[n1, m2] += 1
            cct._C[m1, n1] += 1
        if n2 >= 0:
            cct._B[n2, m2] -= 1
            cct._C[m1, n2] -= 1
        if n3 >= 0:
            cct._B[n3, m1] += 1
            cct._C[m2, n3] += 1
        if n4 >= 0:
            cct._B[n4, m1] -= 1
            cct._C[m2, n4] -= 1                        

        cct._D[m1, m1] += Z2
        cct._D[m2, m2] -= Z1


class CCVS(DependentSource):
    """CCVS"""

    need_branch_current = True
    need_control_current = True

    def _stamp(self, cct):
        n1, n2 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] += 1
            cct._C[m, n1] += 1
        if n2 >= 0:
            cct._B[n2, m] -= 1
            cct._C[m, n2] -= 1
        
        mc = cct._branch_index(self.args[0])
        G = ConstantDomainExpression(self.args[1]).expr
        cct._D[m, mc] -= G

    def _kill(self):
        newopts = self.opts.copy()
        newopts.strip_current_labels()
        newopts.strip_labels()

        return self._netmake_O(opts=newopts)                


class I(IndependentSource):

    add_parallel = True    

    def _select(self, kind=None):
        """Select domain kind for component."""
        return self._netmake(args=self.cpt.Isc.netval(kind))

    def _kill(self):
        newopts = self.opts.copy()
        newopts.strip_voltage_labels()
        newopts.strip_labels()

        return self._netmake_O(opts=newopts)                        

    def _stamp(self, cct):

        n1, n2 = self.node_indexes

        I = self.Isc.expr

        if n1 >= 0:
            cct._Is[n1] += I
        if n2 >= 0:
            cct._Is[n2] -= I

    def _ss_model(self):
        return self._netmake(args='%s(t)' % self.relname.lower())

    def _s_model(self, var):
        return self._netmake(args=self.Isc.laplace()(var))

    def _pre_initial_model(self):

        return self._netmake(args=self.cpt.Isc.pre_initial_value())        


class K(Dummy):
    
    def __init__(self, cct, namespace, defname, name, cpt_type, cpt_id, string,
                 opts_string, nodes, keyword, *args):

        self.Lname1 = args[0]
        self.Lname2 = args[1]
        super (K, self).__init__(cct, namespace, defname, name,
                                 cpt_type, cpt_id, string,
                                 opts_string, nodes, keyword, *args)

    def _stamp(self, cct):

        if cct.kind == 'dc':
            return
        
        if cct.kind in ('t', 'time'):
            raise RuntimeError('Should not be evaluating mutual inductance in'
                               ' time domain')

        L1 = self.Lname1
        L2 = self.Lname2
        K = self.cpt.K

        ZL1 = cct.elements[L1].Z.expr
        ZL2 = cct.elements[L2].Z.expr

        if cct.kind in ('s', 'ivp', 'laplace'):
            # FIXME, generalise for other domains...
            ZM = K.expr * sym.sqrt(ZL1 * ZL2 / s**2) * s
        else:
            ZM = K.expr * sym.sqrt(ZL1 * ZL2)
            
        m1 = cct._branch_index(L1)
        m2 = cct._branch_index(L2)

        cct._D[m1, m2] += -ZM
        cct._D[m2, m1] += -ZM


class L(RLC):
    
    need_branch_current = True
    reactive = True
    add_series = True    

    def _r_model(self):

        dummy_node = self.dummy_node()
        opts = self.opts.copy()        

        # Use Thevenin model.  This will require the current through
        # the voltage source to be explicitly computed.
        
        Req = 'R%seq' % self.name
        Veq = 'V%seq' % self.name        

        opts.strip_voltage_labels()
        rnet = self._netmake_variant('R', suffix='eq',
                                     nodes=(self.relnodes[0], dummy_node),
                                     args=Req, opts=opts)

        opts.strip_current_labels()                
        vnet = self._netmake_variant('V', suffix='eq',
                                     nodes=(dummy_node, self.relnodes[1]),
                                     args=('dc', Veq), opts=opts)

        # TODO: the voltage labels should be added across an
        # open-circuit object.

        return rnet + '\n' + vnet

    @property
    def I0(self):
        """Initial current (for capacitors only)."""

        if self.cct.kind == 'ivp' and self.cpt.has_ic:
            return current(self.cpt.i0 / s)
        return current(0)
    
    @property
    def L(self):
        return self.cpt.L
    
    def _kill(self):
        """Kill implicit sources due to initial conditions."""
        return self.netmake(args=self.args[0])        

    def _initialize(self, ic):
        """Change initial condition to ic."""
        return self._netmake(args=(self.args[0], ic))

    def _stamp(self, cct):

        # This formulation adds the inductor current to the unknowns

        n1, n2 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] = 1
            cct._C[m, n1] = 1
        if n2 >= 0:
            cct._B[n2, m] = -1
            cct._C[m, n2] = -1

        if cct.kind == 'dc':
            Z = 0
        else:
            Z = self.Z.expr

        cct._D[m, m] += -Z

        if cct.kind == 'ivp' and self.cpt.has_ic:
            V = self.Voc.expr            
            cct._Es[m] += V

    def _ss_model(self):
        # Perhaps mangle name to ensure it does not conflict
        # with another current source?
        return self._netmake_variant('I_', args='-i_%s(t)' % self.relname)
    
    def _pre_initial_model(self):

        return self._netmake_variant('I', args=self.cpt.i0)

    
class O(Dummy):
    """Open circuit"""

    def _stamp(self, cct):
        pass

    @property
    def I(self):
        return SuperpositionCurrent(0)

    @property
    def i(self):
        return SuperpositionCurrent(0)(t)
    

class P(O):
    """Port"""
    pass


class R(RC):

    add_series = True    
    
    def _r_model(self):    
        return self._copy()

class NR(R):

    add_series = True    
    
    def _r_model(self):    
        return self._copy()    


class RV(RC):

    # TODO.  Can simulate as series resistors (1 - alpha) R and alpha R. 
    pass


class SPpp(Dummy):

    need_branch_current = True

    def _stamp(self, cct):
        n1, n2, n3 = self.node_indexes
        m = self.branch_index

        if n3 >= 0:
            cct._B[n3, m] += 1
            cct._C[m, n3] += 1
        
        if n1 >= 0:
            cct._C[m, n1] -= 1
        if n2 >= 0:
            cct._C[m, n2] -= 1


class SPpm(Dummy):

    need_branch_current = True

    def _stamp(self, cct):
        n1, n2, n3 = self.node_indexes
        m = self.branch_index

        if n3 >= 0:
            cct._B[n3, m] += 1
            cct._C[m, n3] += 1
        
        if n1 >= 0:
            cct._C[m, n1] -= 1
        if n2 >= 0:
            cct._C[m, n2] += 1

            
class SPppp(Dummy):

    need_branch_current = True

    def _stamp(self, cct):
        n1, n2, n3, n4 = self.node_indexes
        m = self.branch_index

        if n3 >= 0:
            cct._B[n3, m] += 1
            cct._C[m, n3] += 1
        
        if n1 >= 0:
            cct._C[m, n1] -= 1
        if n2 >= 0:
            cct._C[m, n2] -= 1
        if n4 >= 0:
            cct._C[m, n4] -= 1

            
class SPpmm(Dummy):

    need_branch_current = True

    def _stamp(self, cct):
        n1, n2, n3, n4 = self.node_indexes
        m = self.branch_index

        if n3 >= 0:
            cct._B[n3, m] += 1
            cct._C[m, n3] += 1
        
        if n1 >= 0:
            cct._C[m, n1] -= 1
        if n2 >= 0:
            cct._C[m, n2] += 1
        if n4 >= 0:
            cct._C[m, n4] += 1


class SPppm(Dummy):

    need_branch_current = True

    def _stamp(self, cct):
        n1, n2, n3, n4 = self.node_indexes
        m = self.branch_index

        if n3 >= 0:
            cct._B[n3, m] += 1
            cct._C[m, n3] += 1
        
        if n1 >= 0:
            cct._C[m, n1] -= 1
        if n2 >= 0:
            cct._C[m, n2] -= 1
        if n4 >= 0:
            cct._C[m, n4] += 1


class TF(Cpt):
    """Transformer"""    

    need_branch_current = True

    def _stamp(self, cct):

        n1, n2, n3, n4 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] += 1
            cct._C[m, n1] += 1
        if n2 >= 0:
            cct._B[n2, m] -= 1
            cct._C[m, n2] -= 1
        
        T = self.cpt.alpha.expr

        if n3 >= 0:
            cct._B[n3, m] -= T
            cct._C[m, n3] -= T
        if n4 >= 0:
            cct._B[n4, m] += T
            cct._C[m, n4] += T


class TFtap(Cpt):
    """Tapped transformer"""    

    def _stamp(self, cct):
        raise NotImplementedError('Cannot analyse tapped transformer %s' % self)


class TL(Ignored):
    """Transmission line"""

    equipotential_nodes = (('out1', 'in1'),
                           ('out2', 'in2'))
    

class Cable(Ignored):
    """Cable"""

    equipotential_nodes = (('in+', 'out+'), ('in-', 'out-'), ('in', 'out'),
                           ('ignd', 'ognd', 'b', 't'), ('mid', 'out'))


class TP(Misc):
    """Two port"""

    # TODO
    pass


class TR(Dummy):
    """Transfer function.  This is equivalent to a VCVS with the input and
    output referenced to node 0."""

    need_branch_current = True

    def _stamp(self, cct):
        n1, n2 = self.node_indexes
        m = self.branch_index

        if n2 >= 0:
            cct._B[n2, m] += 1
            cct._C[m, n2] += 1
        
        A = ConstantDomainExpression(self.args[0]).expr
        
        if n1 >= 0:
            cct._C[m, n1] -= A


class V(IndependentSource):

    need_branch_current = True
    flip_branch_current = True
    add_series = True

    def _select(self, kind=None):
        """Select domain kind for component."""
        return self._netmake(args=self.cpt.Voc.netval(kind))

    def _kill(self):
        newopts = self.opts.copy()
        newopts.strip_current_labels()
        newopts.strip_labels()

        return self._netmake_W(opts=newopts)

    def _stamp(self, cct):

        n1, n2 = self.node_indexes
        m = self.branch_index

        if n1 >= 0:
            cct._B[n1, m] += 1
            cct._C[m, n1] += 1
        if n2 >= 0:
            cct._B[n2, m] -= 1
            cct._C[m, n2] -= 1

        V = self.Voc.expr
        cct._Es[m] += V

    def _ss_model(self):
        return self._netmake(args='%s(t)' % self.relname.lower())

    def _s_model(self, var):
        return self._netmake(args=self.cpt.Voc.laplace()(var))

    def _pre_initial_model(self):

        return self._netmake(args=self.cpt.Voc.pre_initial_value())


class W(Dummy):
    """Wire"""

    def _stamp(self, cct):
        pass

    @property
    def I(self):
        raise ValueError('Cannot determine current through wire, use a 0 V voltage source')

    @property
    def i(self):
        raise ValueError('Cannot determine current through wire, use a 0 V voltage source')    
    

class XT(Misc):
    """Crystal"""

    reactive = True    


class Y(RC):
    """Admittance"""

    reactive = True
    add_parallel = True


class Z(RC):
    """Impedance"""

    reactive = True
    add_series = True    


classes = {}

def defcpt(name, base, docstring):
    
    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    classes[name] = newclass


def make(classname, parent, name, cpt_type, cpt_id,
         string, opts_string, nodes, *args):

    # Create instance of component object
    newclass = classes[classname]

    cpt = newclass(parent, name, cpt_type, cpt_id, string, opts_string, 
                   nodes, *args)

    return cpt


# Dynamically create classes.
defcpt('A', Misc, 'Annotation')
defcpt('ADC', Misc, 'ADC')
defcpt('AM', W, 'Ammeter')

defcpt('BAT', V, 'Battery')

defcpt('D', NonLinear, 'Diode')
defcpt('DAC', Misc, 'DAC')
defcpt('Dled', 'D', 'LED')
defcpt('Dphoto', 'D', 'Photo diode')
defcpt('Dschottky', 'D', 'Schottky diode')
defcpt('Dtunnel', 'D', 'Tunnel diode')
defcpt('Dzener', 'D', 'Zener diode')

defcpt('E', VCVS, 'VCVS')
defcpt('Eopamp', VCVS, 'Opamp')
defcpt('Efdopamp', VCVS, 'Fully differential opamp')
defcpt('Eamp', VCVS, 'Amplifier')

defcpt('F', CCCS, 'CCCS')
defcpt('FS', Misc, 'Fuse')

defcpt('G', VCCS, 'VCCS')

defcpt('H', CCVS, 'CCVS')

defcpt('sI', I, 's-domain current source')
defcpt('Isin', I, 'Sinusoidal current source')
defcpt('Idc', I, 'DC current source')
defcpt('Istep', I, 'Step current source')
defcpt('Iac', I, 'AC current source')
defcpt('Inoise', I, 'Noise current source')

defcpt('J', NonLinear, 'N JFET transistor')
defcpt('Jnjf', 'J', 'N JFET transistor')
defcpt('Jpjf', 'J', 'P JFET transistor')

defcpt('M', NonLinear, 'N MOSJFET transistor')
defcpt('Mnmos', 'M', 'N channel MOSJFET transistor')
defcpt('Mpmos', 'M', 'P channel MOSJFET transistor')
defcpt('MISC', Misc, 'Generic circuitikz bipole')
defcpt('MT', Misc, 'Motor')
defcpt('MX', Misc, 'Mixer')

defcpt('Q', NonLinear, 'NPN transistor')
defcpt('Qpnp', 'Q', 'PNP transistor')
defcpt('Qnpn', 'Q', 'NPN transistor')

defcpt('Sbox', Misc, 'Box shape')
defcpt('Scircle', Misc, 'Circle shape')
defcpt('Sellipse', Misc, 'Ellipse shape')
defcpt('Striangle', Misc, 'Triangle shape')
defcpt('SW', TimeVarying, 'Switch')
defcpt('SWno', 'SW', 'Normally open switch')
defcpt('SWnc', 'SW', 'Normally closed switch')
defcpt('SWpush', 'SW', 'Pushbutton switch')
defcpt('SWspdt', 'SW', 'SPDT switch')

defcpt('TFcore', TF, 'Transformer with core')
defcpt('TFtapcore', TFtap, 'Transformer with core')

defcpt('Ubuffer', Logic, 'Buffer')
defcpt('Upbuffer', Logic, 'Buffer with power supplies')
defcpt('Uinverter', Logic, 'Inverter')
defcpt('Upinverter', Logic, 'Inverter with power supplies')
defcpt('Uinamp', Misc, 'Instrumentation amplifier')
defcpt('Uisoamp', Misc, 'Isolated amplifier')
defcpt('Udiffamp', Misc, 'Differential amplifier')
defcpt('Udiffdriver', Misc, 'Differential driver')
defcpt('Ufdopamp', Misc, 'Fully differential opamp ')
defcpt('Uopamp', Misc, 'Opamp ')
defcpt('Uregulator', Misc, 'Regulator ')
defcpt('Uadc', Misc, 'ADC')
defcpt('Udac', Misc, 'DAC')
defcpt('Ubox', Misc, 'Box')
defcpt('Ucircle', Misc, 'Circle')
defcpt('Ubox4', Misc, 'Box')
defcpt('Ubox12', Misc, 'Box')
defcpt('Ucircle4', Misc, 'Circle')
defcpt('Uchip1313', Logic, 'General purpose chip')
defcpt('Uchip2121', Logic, 'General purpose chip')
defcpt('Uchip2222', Logic, 'General purpose chip')
defcpt('Uchip3131', Logic, 'General purpose chip')
defcpt('Uchip3333', Logic, 'General purpose chip')
defcpt('Uchip4141', Logic, 'General purpose chip')
defcpt('Uchip4444', Logic, 'General purpose chip')
defcpt('Uchip8181', Logic, 'General purpose chip')
defcpt('Uchip8888', Logic, 'General purpose chip')
defcpt('Umux21', Logic, '2-1 multiplexer')
defcpt('Umux41', Logic, '4-1 multiplexer')
defcpt('Umux42', Logic, '4-2 multiplexer')
defcpt('Udff', Misc, 'D flip-flop')
defcpt('Ujkff', Misc, 'JK flip-flop')
defcpt('Urslatch', Misc, 'RS latch')

defcpt('sV', V, 's-domain voltage source')
defcpt('Vsin', V, 'Sinusoidal voltage source')
defcpt('Vdc', V, 'DC voltage source')
defcpt('Vstep', V, 'Step voltage source')
defcpt('Vac', V, 'AC voltage source')
defcpt('Vnoise', V, 'Noise voltage source')
defcpt('VM', O, 'Voltmeter')

# Let's choose mechanical analogue II (the impedance analogue) where
# force is equivalent to voltage and velocity is equivalent to
# current.  With this analogy parallel and serial have to be switched.
defcpt('m', L, 'Mass')
defcpt('k', C, 'Spring')
defcpt('r', R, 'Damper')

# Append classes defined in this module but not imported.
clsmembers = inspect.getmembers(module, lambda member: inspect.isclass(member) and member.__module__ == __name__)
for name, cls in clsmembers:
    classes[name] = cls

