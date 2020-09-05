import numpy as np
from .oneport import Ser, Par
from .netlisthelper import NetlistHelper

class LadderMaker(NetlistHelper):

    def __init__(self, net):

        self.net = net
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

        # Have Par or Ser.  Only one arg can have a depth more than 1.
        depths = net._depths
        num = sum([depth > 0 for depth in depths])
        if num > 1:
            raise ValueError('Cannot draw as ladder network')

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
                self._add('W %s %s; right\n' % (n2, n2p))
                self._add('W %s %s; down\n' % (n1p, n2p))
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
            self._add('W %s %s; right\n' % (n2, n2p))            

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

        # Add wires at start if Par
        if isinstance(net, Par):
            n1p = self._node
            n2p = self._node
            self._add('W %s %s; right=0.5\nW %s %s; right=0.5\n' % (n2, n2p, n1, n1p))
            n1, n2 = n1p, n2p

        depths = net._depths
        if len(depths) != 2:
            # FIXME  Might have Par(N1, N2, N3)
            raise ValueError('Cannot draw as ladder network')

        self.first_cpt = net.args[self._min_depth(depths)]

        self._section_make(net, n1, n2)
        return self.s
    
