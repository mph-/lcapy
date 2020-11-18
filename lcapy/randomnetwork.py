"""This module contains the RandomNetwork class for creating random
networks.

Copyright 2020 Michael Hayes, UCECE"""

from .oneport import Vdc, Idc, Vstep, Istep, Vac, Iac, R, L, C
import random


class RandomNetworkMaker(object):

    def __init__(self, NR=3, NL=0, NC=0, NV=1, NI=0,
                 Nparallel=None, kind='transient'):
        
        self.NR = NR
        self.NL = NL
        self.NC = NC
        self.NV = NV
        self.NI = NI
        self.Nparallel = Nparallel
        self.kind = kind

        cpts = []
        cpts.extend(self._add_cpts(R, self.NR))
        cpts.extend(self._add_cpts(L, self.NL))
        cpts.extend(self._add_cpts(C, self.NC))

        if self.kind == 'transient':
            cpts.extend(self._add_cpts(Vstep, self.NV))
            cpts.extend(self._add_cpts(Istep, self.NI))
        elif self.kind == 'dc':
            cpts.extend(self._add_cpts(Vdc, self.NV))
            cpts.extend(self._add_cpts(Idc, self.NI))
        elif self.kind == 'ac':
            cpts.extend(self._add_cpts(Vac, self.NV))
            cpts.extend(self._add_cpts(Iac, self.NI))            
        else:
            raise ValueError('Unknown circuit kind %s' % self.kind)

        self.cpts = cpts
        self.Nparallel = Nparallel

        self.connections = self._choose_connections()

    def _choose_connections(self):
        Nconnections = len(self.cpts) - 1
        
        if self.Nparallel is None:
            connections = random.choices((True, False), k=Nconnections)
        else:
            Nparallel = min(self.Nparallel, Nconnections)
            Nseries = Nconnections - Nparallel
            connections = [True] * Nparallel + [False] * Nseries

        return connections

    def __call__(self):

        return self.make()

    def _make1(self):

        cpts = random.sample(self.cpts, len(self.cpts))
        connections = random.sample(self.connections, len(self.connections))        

        net = cpts[0]
        for cpt, connection in zip(cpts[1:], connections):
            if connection:
                if cpt.voltage_source and net.has_parallel_V:
                    return None
                net = net.parallel(cpt)
            else:
                if cpt.current_source and net.has_series_I:
                    return None
                net = net.series(cpt)            
        return net

    def make(self):

        while True:
            net = self._make1()
            if net is not None:
                return net
            self.connections = self._choose_connections()
            
    def _add_cpts(self, cpt, num):

        cpts = []
        for m in range(num):
            cpts.append(cpt('%s%d' % (cpt.__name__[0:1], m + 1)))        
        return cpts

    
def random_network(NR=3, NL=0, NC=0, NV=1, NI=0, Nparallel=None,
                   kind='transient'):

    return RandomNetworkMaker(NR=NR, NL=NL, NC=NC, NV=NV, NI=NI,
                              Nparallel=Nparallel, kind=kind).make()
