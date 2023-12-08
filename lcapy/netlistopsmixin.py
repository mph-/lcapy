"""This module provides the NetlistOpsMixin class.

Copyright 2023 Michael Hayes, UCECE

"""

# TODO: Optimise oneport in a similar manner to ladder.  This will avoid
# solving for all nodal voltages.  Voc, ISc, impedance, and admittance
# can then be sped up using the oneport.

from .admittance import admittance
from .current import current, current_sign
from .impedance import impedance
from .matrix import Matrix
from .statespace import StateSpace
from .symbols import s
from .transfer import transfer
from warnings import warn


class NetlistOpsMixin:

    def admittance(self, Np, Nm=None):
        """Return driving-point Laplace-domain admittance between nodes
        Np and Nm with independent sources killed and initial
        conditions ignored."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.kill()
        new._add_ground(Nm)
        test = new._add_test_voltage_source(Np, Nm)
        If = current_sign(-new[test].I, True)
        return admittance(If.laplace().sympy)

    def conductance(self, Np, Nm=None):
        """Return conductance (inverse resistance) between nodes Np and Nm
        with independent sources killed.  The result is in the AC (omega)
        domain.    See also resistance, reactance, susceptance."""
        return self.impedance(Np, Nm).G

    def current_gain(self, N1p, N1m, N2p=None, N2m=None):
        """Create Laplace-domain current transfer function I2(s) / I1(s) where:
        I1 is the test current applied between N1p and N1m
        I2 is the measured short-circuit current flowing from N2m to N2p

        Note, the currents are considered to be flowing into the
        positive nodes as is the convention with two-ports.  Thus the
        input and output currents have opposite directions and so a
        piece of wire has a current gain of -1.

        Note, independent sources are killed and initial conditions
        are ignored.

        Alternative forms are:
            current_gain(N1p, N1m, N2p, N2m)
            current_gain(cpt1, cpt2)
            current_gain((N1p, N1m), cpt2)
            current_gain(cpt1, (N2p, N2m))

        """

        try:
            # Finding the transfer function for a ladder network
            # is faster but is not necessary
            ladder = self._ladder(N1p, N1m, N2p, N2m)
            if ladder is not None:
                return ladder.current_gain
        except:
            pass

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'current_gain')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_current_source(N1p, N1m)
        H = transfer(-new.Isc(N2p, N2m).laplace())
        H.causal = True
        return H

    def impedance(self, Np, Nm=None):
        """Return driving-point Laplace-domain impedance between nodes
        Np and Nm with independent sources killed and initial
        conditions ignored."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.apply_test_current_source(Np, Nm)
        Vf = new.Voc(Np, Nm)
        return impedance(Vf.laplace().sympy)

    def Isc(self, Np, Nm=None, **kwargs):
        """Return short-circuit transform-domain current between nodes Np and
        Nm."""

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)

        new = self.copy()
        if new.is_causal:
            new.add('Vshort_ %s %s step 0' % (Np, Nm))
        else:
            new.add('Vshort_ %s %s 0' % (Np, Nm))

        Isc = current_sign(new.get_I('Vshort_', **kwargs), True)

        new.remove('Vshort_')

        return Isc

    def isc(self, Np, Nm=None):
        """Return short-circuit t-domain current between nodes Np and Nm."""

        return self.Isc(Np, Nm).time()

    def _ladder(self, N1p, N1m, N2p=None, N2m=None):

        # TODO: determine point at which it is more efficient to
        # find ladder network
        if len(self.elements) < 6:
            return None
        return self.kill().ladder(N1p, N1m, N2p, N2m)

    def ladder(self, N1p, N1m, N2p=None, N2m=None):
        """Return two-port unbalanced ladder network or `None` if the netlist
        does not have a ladder topology between the specified nodes.

        The input port is defined by the nodes `N1p` and `N1m`.
        The output port is defined by the nodes `N2p` and `N2m`.

        The nodes `N1p` and `N1m` must be the same.

        Alternative forms are:
            ladder(N1p, N1m, N2p, N2m)
            ladder(cpt1, cpt2)
            ladder((N1p, N1m), cpt2)
            ladder(cpt1, (N2p, N2m))

        This method can be used to generate a transfer function, for example:

        `cct.ladder(1, 0, 10, 0).voltage_gain`

        This can be much faster than using `cct.voltage_gain(1, 0, 10, 0)`
        for circuits with a large ladder topology, since the latter method needs
        to determine all the node voltages.

        """

        from lcapy.laddernetworkmaker import LadderNetworkMaker

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'ladder')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        lm = LadderNetworkMaker(self)

        tp = lm.make(N1p, N1m, N2p, N2m)
        return tp

    def norton(self, Np, Nm=None):
        """Return Laplace-domain Norton model between nodes Np and Nm.

        If Np is a component name, create model using its component nodes."""

        from .oneport import I, Y

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)
        Isc = self.Isc(Np, Nm)
        Ysc = self.admittance(Np, Nm)

        return (I(Isc) | Y(Ysc)).simplify()

    def oneport(self, Np, Nm=None):
        """Return oneport object between nodes Np and Nm.  This might be a
        Thevenin network, a Norton network, or a single component.

        If Np is a component name, create model using its component nodes."""

        try:
            return self.norton(Np, Nm)
        except:
            return self.thevenin(Np, Nm)

    def reactance(self, Np, Nm=None):
        """Return reactance between nodes Np and Nm with independent
        sources killed.  The result is in the AC (omega) domain.
        See also conductance, resistance, susceptance."""
        return self.impedance(Np, Nm).X

    def resistance(self, Np, Nm=None):
        """Return resistance between nodes Np and Nm with independent
        sources killed.  The result is in the AC (omega) domain.
        See also conductance, reactance, susceptance."""
        return self.impedance(Np, Nm).R

    def state_space(self, node_voltages=None, branch_currents=None):
        """Generate state-space representation.

        `node_voltages` is a list of node names to use as voltage outputs.
        If `None` use all the unique node names.

        `branch_currents` is a list of component names to use as
        current outputs.  If `None` use all the components.

        Here's an example:
        `cct = Circuit('cct.sch')
        ss = cct.state_space(node_voltages=['1', '3'], branch_currents=['L1', 'L2'])`
        """

        ss = StateSpace.from_circuit(self, node_voltages, branch_currents)
        return ss

    def susceptance(self, Np, Nm=None):
        """Return susceptance (inverse reactance) between nodes Np and Nm with
        independent sources killed.  The result is in the AC (omega)
        domain.  See also conductance, reactance, resistance."""
        return self.impedance(Np, Nm).B

    def thevenin(self, Np, Nm=None):
        """Return Laplace-domain Thevenin oneport model between nodes Np and Nm.

        If Np is a component name, create model using its component nodes."""

        from .oneport import V, Z

        Np, Nm = self._parse_node_args2(Np, Nm)
        Np, Nm = self._check_nodes(Np, Nm)
        Voc = self.Voc(Np, Nm)
        Zoc = self.impedance(Np, Nm)

        return (V(Voc) + Z(Zoc)).simplify()

    def transfer(self, N1p, N1m, N2p=None, N2m=None):
        """Create Laplace-domain voltage transfer function V2(s) / V1(s) where:
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        Note, independent sources are killed and initial conditions
        are ignored.

        Alternative forms are:
            transfer(N1p, N1m, N2p, N2m)
            transfer(cpt1, cpt2)
            transfer((N1p, N1m), cpt2)
            transfer(cpt1, (N2p, N2m))
        """

        try:
            # Finding the transfer function for a ladder network
            # is faster but is not necessary
            ladder = self._ladder(N1p, N1m, N2p, N2m)
            if ladder is not None:
                return ladder.voltage_gain
        except:
            pass

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'transfer')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_voltage_source(N1p, N1m)
        V2 = new.Voc(N2p, N2m)
        H = transfer(V2.laplace())
        H.causal = True
        return H

    def twoport(self, N1p, N1m, N2p=None, N2m=None, model='B'):
        """Create Laplace-domain twoport model for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        `model` is `A, `B`, `G`, `H`, `Y`, or `Z`.

        Alternative forms are:
            twoport(N1p, N1m, N2p, N2m)
            twoport(cpt1, cpt2)
            twoport((N1p, N1m), cpt2)
            twoport(cpt1, (N2p, N2m))
        """

        from .twoport import TwoPortAModel, TwoPortBModel, TwoPortGModel
        from .twoport import TwoPortHModel, TwoPortYModel, TwoPortZModel

        N1p, N1m, N2p, N2m = self._parse_node_args4(
            N1p, N1m, N2p, N2m, 'twoport')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        # TODO, generalise for not just Laplace-domain.

        new = self.copy()
        new._add_ground(N1m)

        if model == 'A':
            V1a = new.Voc(N1p, N1m, nowarn=True)(s)
            I1a = new.Isc(N1p, N1m, nowarn=True)(s)
            A = new.Aparams(N1p, N1m, N2p, N2m)
            return TwoPortAModel(A, V1a=V1a, I1a=I1a)
        elif model == 'B':
            V2b = new.Voc(N2p, N2m, nowarn=True)(s)
            I2b = new.Isc(N2p, N2m, nowarn=True)(s)
            A = new.Aparams(N1p, N1m, N2p, N2m)
            return TwoPortBModel(A.Bparams, V2b=V2b, I2b=I2b)
        elif model == 'Z':
            V1 = new.Voc(N1p, N1m, nowarn=True)(s)
            V2 = new.Voc(N2p, N2m, nowarn=True)(s)
            Z = new.Zparams(N1p, N1m, N2p, N2m)
            return TwoPortZModel(Z, V1z=V1, V2z=V2)
        elif model == 'Y':
            I1 = new.Isc(N1p, N1m, nowarn=True)(s)
            I2 = new.Isc(N2p, N2m, nowarn=True)(s)
            Z = new.Zparams(N1p, N1m, N2p, N2m)
            return TwoPortYModel(Z.Yparams, I1y=I1, I2y=I2)
        elif model == 'G':
            I1 = new.Isc(N1p, N1m, nowarn=True)(s)
            V2 = new.Voc(N2p, N2m, nowarn=True)(s)
            Z = new.Zparams(N1p, N1m, N2p, N2m)
            return TwoPortGModel(Z.Gparams, I1g=I1, V2g=V2)
        elif model == 'H':
            V1 = new.Voc(N1p, N1m, nowarn=True)(s)
            I2 = new.Isc(N2p, N2m, nowarn=True)(s)
            Z = new.Zparams(N1p, N1m, N2p, N2m)
            return TwoPortHModel(Z.Hparams, V1h=V1, I2h=I2)
        else:
            raise ValueError('Model %s unknown, must be B, H, Y, or Z' % model)

    def transadmittance(self, N1p, N1m, N2p=None, N2m=None):
        """Create Laplace-domain transadmittance (transfer admittance) function
        I2(s) / V1(s) where:
          V1 is the test voltage applied between N1p and N1m
          I2 is the measured short-circuit current flowing from N2m to N2p.

        Note, I2 is considered to be flowing into the positive node as
        is the convention with two-ports.  Thus the transadmittance of
        a series resistor with resistance R is -1 / R.

        Note, independent sources are killed and initial conditions
        are ignored.

        Alternative forms are:
            transadmittance(N1p, N1m, N2p, N2m)
            transadmittance(cpt1, cpt2)
            transadmittance((N1p, N1m), cpt2)
            transadmittance(cpt1, (N2p, N2m))

        """

        try:
            # Finding the transfer function for a ladder network
            # is faster but is not necessary
            ladder = self._ladder(N1p, N1m, N2p, N2m)
            if ladder is not None:
                return ladder.transadmittance
        except:
            pass

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'transadmittance')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_voltage_source(N1p, N1m)
        H = admittance(-new.Isc(N2p, N2m).laplace())
        H.causal = True
        return H

    def transimpedance(self, N1p, N1m, N2p=None, N2m=None):
        """Create Laplace-domain transimpedance (transfer impedance) function
        V2(s) / I1(s) where:
          I1 is the test current applied between N1p and N1m
          V2 is the measured open-circuit voltage between N2p and N2m.

        Note, I1 is considered to be flowing into the positive node as
        is the convention with two-ports.

        Note, independent sources are killed and initial conditions
        are ignored.

        Alternative forms are:
            transimpedance(N1p, N1m, N2p, N2m)
            transimpedance(cpt1, cpt2)
            transimpedance((N1p, N1m), cpt2)
            transimpedance(cpt1, (N2p, N2m))
        """

        try:
            # Finding the transfer function for a ladder network
            # is faster but is not necessary
            ladder = self._ladder(N1p, N1m, N2p, N2m)
            if ladder is not None:
                return ladder.transimpedance
        except:
            pass

        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'transadmittance')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_current_source(N1p, N1m)
        H = impedance(new.Voc(N2p, N2m).laplace())
        H.causal = True
        return H

    def Voc(self, Np, Nm=None, **kwargs):
        """Return open-circuit transform-domain voltage between nodes Np and
        Nm."""

        return self._get_Vd(Np, Nm, **kwargs)

    def voc(self, Np, Nm=None):
        """Return open-circuit t-domain voltage between nodes Np and Nm."""

        return self.Voc(Np, Nm).time()

    def voltage_gain(self, N1p, N1m, N2p=None, N2m=None):
        """Create Laplace-domain voltage transfer function V2(s) / V1(s) where:
        V1 is the test voltage applied between N1p and N1m
        V2 is the measured open-circuit voltage between N2p and N2m

        Note, independent sources are killed and initial conditions
        are ignored.

        Alternative forms are:
            voltage_gain(N1p, N1m, N2p, N2m)
            voltage_gain(cpt1, cpt2)
            voltage_gain((N1p, N1m), cpt2)
            voltage_gain(cpt1, (N2p, N2m))
        """

        try:
            # Finding the transfer function for a ladder network
            # is faster but is not necessary
            ladder = self._ladder(N1p, N1m, N2p, N2m)
            if ladder is not None:
                return ladder.voltage_gain
        except:
            pass

        # This is the same as transfer.
        N1p, N1m, N2p, N2m = self._parse_node_args4(N1p, N1m, N2p, N2m,
                                                    'voltage_gain')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)

        new = self.apply_test_voltage_source(N1p, N1m)
        V2 = new.Voc(N2p, N2m)
        H = transfer(V2.laplace())
        H.causal = True
        return H

    def Aparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create A-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Bparams, Gparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """

        from .twoport import AMatrix

        N1p, N1m, N2p, N2m = self._parse_node_args4(
            N1p, N1m, N2p, N2m, 'Aparams')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)
        new = self.kill()
        new._add_ground(N1m)

        try:
            test = new._add_test_voltage_source(N1p, N1m)

            # A11 = V1 / V2 with I2 = 0
            # Apply V1 and measure V2 with port 2 open-circuit
            A11 = new.Voc(N1p, N1m)(s) / new.Voc(N2p, N2m)(s)

            # A12 = V1 / I2 with V2 = 0
            # Apply V1 and measure -I2 with port 2 short-circuit
            A12 = new.Voc(N1p, N1m)(s) / new.Isc(N2p, N2m)(s)

            new.remove(test)

            test = new._add_test_current_source(N1p, N1m)

            # A21 = I1 / V2 with I2 = 0
            # Apply I1 and measure V2 with port 2 open-circuit
            try:
                A21 = current(0 * s + 1) / new.Voc(N2p, N2m)(s)
            except ValueError:
                # It is likely there is an open-circuit.
                new2 = new.copy()
                new2.add('W %s %s' % (N2p, N2m))
                A21 = -new2[test].I(s) / new2.Voc(N2p, N2m)(s)
                A21 = 0

            # A22 = I1 / I2 with V2 = 0
            # Apply I1 and measure -I2 with port 2 short-circuit
            A22 = current(0 * s + 1) / new.Isc(N2p, N2m)(s)

            new.remove(test)
            A = AMatrix(((A11, A12), (A21, A22)))
            return A

        except ValueError:
            warn('Cannot create A matrix directly; trying via Z matrix')
            Z = self.Zparams(N1p, N1m, N2p, N2m)
            return Z.Aparams

    def Bparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create B-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Gparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Bparams

    def Gparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create G-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Hparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Gparams

    def Hparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create H-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Sparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Hparams

    def Sparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create S-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Hparams, Tparams, Yparams, and Zparams.
        """
        return self.Aparams(N1p, N1m, N2p, N2m).Sparams

    def Tparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create T-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Hparams, Sparams, Yparams, and Zparams.
        """
        return self.Tparams(N1p, N1m, N2p, N2m).Hparams

    def Yparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create Y-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Hparams, Sparams, Tparams, and Zparams.
        """
        return self.Zparams(N1p, N1m, N2p, N2m).Yparams

    def Zparams(self, N1p, N1m, N2p=None, N2m=None):
        """Create Z-parameters for two-port defined by nodes N1p, N1m, N2p, and N2m, where:
        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]

        See also twoport, Aparams, Bparams, Gparams, Hparams, Sparams, Tparams, and Yparams.
        """
        from .twoport import ZMatrix

        # TODO, generalise to multiports.

        N1p, N1m, N2p, N2m = self._parse_node_args4(
            N1p, N1m, N2p, N2m, 'Zparams')
        N1p, N1m, N2p, N2m = self._check_nodes(N1p, N1m, N2p, N2m)
        new = self.kill()
        new._add_ground(N1m)

        try:
            test = new._add_test_current_source(N1p, N1m)

            # Z11 = V1 / I1 with I2 = 0
            # Apply I1 and measure V1 with port 2 open-circuit
            Z11 = impedance(new.Voc(N1p, N1m)(s))

            # Z21 = V2 / I1 with I2 = 0
            # Apply I1 and measure V2 with port 2 open-circuit
            Z21 = impedance(new.Voc(N2p, N2m)(s))

            new.remove(test)

            test = new._add_test_current_source(N2p, N2m)

            # Z12 = V1 / I2 with I1 = 0
            # Apply I2 and measure V1 with port 1 open-circuit
            Z12 = impedance(new.Voc(N1p, N1m)(s))

            # Z22 = V2 / I2 with I1 = 0
            # Apply I2 and measure V2 with port 1 open-circuit
            Z22 = impedance(new.Voc(N2p, N2m)(s))

            new.remove(test)

            Z = ZMatrix(((Z11, Z12), (Z21, Z22)))
            return Z

        except ValueError as e:
            raise ValueError('Cannot create Z matrix: %s' % e)

    def Yparamsn(self, *nodes):
        """Create Y-parameters for N-port defined by list of node-pairs.

        See also Yparams for a two port.

        """

        nodes = self._check_nodes(*nodes)
        if len(nodes) % 2 == 1:
            raise ValueError('Need an even number of nodes.')
        ports = []
        for m in range(len(nodes) // 2):
            ports.append((nodes[m * 2], nodes[m * 2 + 1]))

        new = self.kill()
        new._add_ground(nodes[1])

        try:

            Y = Matrix.zeros(len(ports))

            for col in range(len(ports)):

                for row in range(len(ports)):
                    if row == col:
                        new.add('V%d_ %s %s {DiracDelta(t)}' % (
                            row, ports[row][0], ports[row][1]))
                    else:
                        new.add('V%d_ %s %s 0' %
                                (row, ports[row][0], ports[row][1]))

                for row in range(len(ports)):
                    Y[row, col] = admittance(new.elements['V%d_' % row].I(s))

                for row in range(len(ports)):
                    new.remove('V%d_' % row)
            return Y

        except ValueError as e:
            raise ValueError('Cannot create Y matrix: %s' % e)

    def Yparams3(self, N1p, N1m, N2p, N2m, N3p, N3m):
        """Create Y-parameters for three-port defined by nodes N1p, N1m, N2p,
        N2m, N3p, and N3m, where:

        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        I3 is the current flowing into N3p and out of N3m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        V3 is V[N3p] - V[N3m]

        See also Yparams for a two port and Yparamsn for an N-port.

        """

        return self.Yparamsn(N1p, N1m, N2p, N2m, N3p, N3m)

    def Zparamsn(self, *nodes):
        """Create Z-parameters for N-port defined by list of node-pairs.

        See also Zparams for a two port.

        """

        nodes = self._check_nodes(*nodes)
        if len(nodes) % 2 == 1:
            raise ValueError('Need an even number of nodes.')
        ports = []
        for m in range(len(nodes) // 2):
            ports.append((nodes[m * 2], nodes[m * 2 + 1]))

        new = self.kill()
        new._add_ground(nodes[1])

        try:

            Z = Matrix.zeros(len(ports))

            for col in range(len(ports)):
                new.add('I_ %s %s {DiracDelta(t)}' %
                        (ports[col][0], ports[col][1]))

                for row in range(len(ports)):
                    Z[row, col] = impedance(
                        new.Voc(ports[row][0], ports[row][1])(s))

                new.remove('I_')
            return Z

        except ValueError as e:
            raise ValueError('Cannot create Z matrix: %s' % e)

    def Zparams3(self, N1p, N1m, N2p, N2m, N3p, N3m):
        """Create Z-parameters for three-port defined by nodes N1p, N1m, N2p,
        N2m, N3p, and N3m, where:

        I1 is the current flowing into N1p and out of N1m
        I2 is the current flowing into N2p and out of N2m
        I3 is the current flowing into N3p and out of N3m
        V1 is V[N1p] - V[N1m]
        V2 is V[N2p] - V[N2m]
        V3 is V[N3p] - V[N3m]

        See also Zparams for a two port and Zparamsn for an N-port.

        """

        return self.Zparamsn(N1p, N1m, N2p, N2m, N3p, N3m)
