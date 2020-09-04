import numpy as np
from .oneport import Ser, Par
from .netlisthelper import NetlistHelper

class LadderMaker(NetlistHelper):

    def __init__(self, net):

        self.net = net
        self.s = ''

    def _add(self, s):

        if not s.endswith('\n'):
            s = s + '\n'
        
        self.s = self.s + s

    def _net_add(self, net, n1, n2, dir='right'):
        self._add(net._net_make(self, n1, n2, dir))
        
    def _section_make(self, net, n1, n2):

        if not isinstance(net, (Ser, Par)):
            return self._net_add(net, n1, n2, dir='down')

        # Have Par or Ser.  Only one arg can have a depth more than 1.
        depths = net._depths
        num = sum([depth > 0 for depth in depths])
        if num > 1:
            raise ValueError('Cannot draw as ladder network')

        if isinstance(net, Par):
            # TODO handle last Par in AST.   This has a depth 1.
            
            mm = None
            for m, net1 in enumerate(net.args):
                if depths[m] != 0:
                    mm = m
                    continue
                # TODO, for parallel, need extra wires
                self._net_add(net1, n1, n2, dir='down')

            if mm is None:
                return
            return self._section_make(net.args[mm], n1, n2)

        elif isinstance(net, Ser):
            # TODO handle last Ser in AST.   This has a depth 1.

            # TODO handle balanced/unbalanced by halving series
            # impedance and inserting in each branch.
            
            n1p = self._node
            n2p = self._node            
            
            mm = None
            for m, net1 in enumerate(net.args):
                if depths[m] != 0:
                    mm = m
                    continue
                self._net_add(net1, n1, n1p, dir='right')

            self._add('W %s %s; dir=right\n' % (n2, n2p))
                
            if mm is None:
                return
            return self._section_make(net.args[mm], n1p, n2p)

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
    
