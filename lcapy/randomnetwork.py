"""This module contains the RandomNetwork class for creating random
networks.

Copyright 2020--2021 Michael Hayes, UCECE"""

from .oneport import Vdc, Idc, Vstep, Istep, Vac, Iac, R, L, C
import random


class RandomNetworkMaker(object):

    def __init__(self, num_resistors=3, num_inductors=0, num_capacitors=0,
                 num_voltage_sources=1, num_current_sources=0,
                 num_parallel=None, numeric_values=False, kind='transient'):
        """`kind` can be 'transient', 'ac', or 'dc'."""
        
        self.num_resistors = num_resistors
        self.num_inductors = num_inductors
        self.num_capacitors = num_capacitors
        self.num_voltage_sources = num_voltage_sources
        self.num_current_sources = num_current_sources
        self.num_parallel = num_parallel
        self.numeric_values = numeric_values
        self.kind = kind

        cpts = []
        cpts.extend(self._add_cpts(R, self.num_resistors))
        cpts.extend(self._add_cpts(L, self.num_inductors))
        cpts.extend(self._add_cpts(C, self.num_capacitors))

        if self.kind == 'transient':
            cpts.extend(self._add_cpts(Vstep, self.num_voltage_sources))
            cpts.extend(self._add_cpts(Istep, self.num_current_sources))
        elif self.kind == 'dc':
            cpts.extend(self._add_cpts(Vdc, self.num_voltage_sources))
            cpts.extend(self._add_cpts(Idc, self.num_current_sources))
        elif self.kind == 'ac':
            cpts.extend(self._add_cpts(Vac, self.num_voltage_sources))
            cpts.extend(self._add_cpts(Iac, self.num_current_sources))            
        else:
            raise ValueError('Unknown circuit kind %s' % self.kind)

        self.cpts = cpts
        self.num_parallel = num_parallel

        self.connections = self._choose_connections()

    def _choose_connections(self):
        Nconnections = len(self.cpts) - 1
        
        if self.num_parallel is None:
            connections = random.choices((True, False), k=Nconnections)
        else:
            num_parallel = min(self.num_parallel, Nconnections)
            Nseries = Nconnections - num_parallel
            connections = [True] * num_parallel + [False] * Nseries

        return connections

    def __call__(self):

        return self.make()

    def _make1(self):

        cpts = random.sample(self.cpts, len(self.cpts))
        connections = random.sample(self.connections, len(self.connections))

        net = cpts[0]
        for cpt, connection in zip(cpts[1:], connections):
            if connection:
                if cpt.is_voltage_source and net.has_parallel_V:
                    return None
                net = net.parallel(cpt)
            else:
                if cpt.is_current_source and net.has_series_I:
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
            if self.numeric_values:
                value = str(random.randint(1, 10))
            else:
                value = '%s%d' % (cpt.__name__[0:1], m + 1)
            
            cpts.append(cpt(value))
        return cpts

    
def random_network(num_resistors=3, num_inductors=0, num_capacitors=0,
                   num_voltage_sources=1, num_current_sources=0,
                   num_parallel=None, numeric_values=False, kind='transient'):

    return RandomNetworkMaker(num_resistors=num_resistors,
                              num_inductors=num_inductors,
                              num_capacitors=num_capacitors,
                              num_voltage_sources=num_voltage_sources,
                              num_current_sources=num_current_sources,
                              num_parallel=num_parallel,
                              numeric_values=numeric_values,
                              kind=kind).make()
