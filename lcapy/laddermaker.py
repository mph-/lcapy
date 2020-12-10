"""This module provides the LadderMaker class for generating
a netlist with a ladder structure from a network.

Copyright 2020 Michael Hayes, UCECE

"""

import numpy as np
from .oneport import Ser, Par
from .netlisthelper import NetlistHelper

class LadderMaker(NetlistHelper):

    def __init__(self, net, form='ladder', evalf=None):

        self.net = net
        self.form = form
        self.evalf = evalf
        self.s = ''
        self.first_cpt = None
        self.first_series = None
        self.first_parallel = None

    def _add(self, s):

        if not s.endswith('\n'):
            s = s + '\n'
        
        self.s = self.s + s

    def _net_add(self, net, n1, n2, dir='right'):
        self._add(net._net_make(self, n1, n2, dir))

    def _wire_add(self, n1, n2, dir='right'):
        self._add('W %s %s; %s\n' % (n1, n2, dir))

    def _port_add(self, n1, n2, dir='right'):
        self._add('P %s %s; %s\n' % (n1, n2, dir))        

    def _split(self, net, cpt):

        if isinstance(net.args[0], cpt.__class__):
            return net.args[0], net.args[1]
        if isinstance(net.args[1], cpt.__class__):
            return net.args[1], net.args[0]

        # Do not have symmetry so choose...
        print('Missing ladder symmetry')
        return net.args[0], net.args[1]        
        
    def _section_make(self, net, n1, n2):

        # Should have alternating Par, Ser
        
        if not isinstance(net, (Ser, Par)):
            return self._net_add(net, n1, n2, dir='down')

        # Have Par or Ser.
        depths = net._depths

        if isinstance(net, Par):
            if len(net.args) > 2:
                raise ValueError('FIXME for more than two parallel cpts')

            if self.first_parallel is None:
                self.first_parallel = net.args[self._min_depth(depths)]

            cpt, rest = self._split(net, self.first_parallel)
            self._net_add(cpt, n1, n2, dir='down')

            if not isinstance(rest, (Ser, Par)):
                n1p = self._node
                n2p = self._node                                
                self._net_add(rest, n1, n1p, dir='right')
                self._wire_add(n2, n2p, 'right')
                self._wire_add(n1p, n2p, 'down')
                return

            return self._section_make(rest, n1, n2)

        elif isinstance(net, Ser):

            if self.first_series is None:
                self.first_series = net.args[self._min_depth(depths)]
            
            # TODO handle balanced/unbalanced by halving series
            # impedance and inserting in each branch.
            
            n1p = self._node
            n2p = self._node

            cpt, rest = self._split(net, self.first_series)
            self._net_add(cpt, n1, n1p, dir='right')
            self._wire_add(n2, n2p, 'right')

            return self._section_make(rest, n1p, n2p)            

        else:
            raise ValueError('Unhandled component %s' % net)
        
    def _min_depth(self, depths):

        return np.argmin(depths)

    def __call__(self):

        net = self.net
        n2 = self._node           # 0
        n1 = self._node           # 1

        if not isinstance(net, (Ser, Par)):
            return self._net_add(net, n1, n2, dir='down')

        self._port_add(n1, n2, dir='down')        

        # Add wires at start if Par
        if isinstance(net, Par):
            n1p = self._node
            n2p = self._node
            self._wire_add(n1, n1p, 'right=0.5')            
            self._wire_add(n2, n2p, 'right=0.5')
            n1, n2 = n1p, n2p

        depths = net._depths
        if len(depths) != 2:
            # FIXME  Might have Par(N1, N2, N3)
            raise ValueError('Cannot draw as ladder network')

        self.first_cpt = net.args[self._min_depth(depths)]

        self._section_make(net, n1, n2)
        return self.s
    
