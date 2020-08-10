"""This module contains the RandomNetwork class for creating random
networks.

Copyright 2020 Michael Hayes, UCECE"""

from .oneport import V, I, Vstep, Istep, Vac, Iac, R, L, C
import random


class RandomNetworkMaker(object):

    def __init__(self, NR=3, NL=0, NC=0, NV=1, NI=0, kind='transient'):
        
        self.NR = NR
        self.NL = NL
        self.NC = NC
        self.NV = NV
        self.NI = NI
        self.kind = kind

    def __call__(self):

        return self.make()

    def make(self):

        cpts = []
        cpts.extend(self._add_cpts(R, self.NR))
        cpts.extend(self._add_cpts(L, self.NL))
        cpts.extend(self._add_cpts(C, self.NC))

        if self.kind == 'transient':
            cpts.extend(self._add_cpts(Vstep, self.NV))
            cpts.extend(self._add_cpts(Istep, self.NI))
        elif kind == 'dc':
            cpts.extend(self._add_cpts(Vstep, self.NV))
            cpts.extend(self._add_cpts(Istep, self.NI))
        elif kind == 'ac':
            cpts.extend(self._add_cpts(Vac, self.NV))
            cpts.extend(self._add_cpts(Iac, self.NI))            
        else:
            raise ValueError('Unknown circuit kind %s' % self.kind)

        cpts = random.sample(cpts, len(cpts))

        net = cpts[0]
        for cpt in cpts[1:]:
            if random.choice((True, False)):
                net = net.parallel(cpt)
            else:
                net = net.series(cpt)            
        return net        
        
    def _add_cpts(self, cpt, num):

        cpts = []
        for m in range(num):
            cpts.append(cpt('%s%d' % (cpt.__name__[0:1], m + 1)))        
        return cpts

    
def random_network(NR=3, NL=0, NC=0, NV=1, NI=0, kind='transient'):

    return RandomNetworkMaker(NR, NL, NV, NV, NI, kind).make()


