"""
This module collects information about a circuit.

Copyright 2022 Michael Hayes, UCECE

"""


class Analysis:

    def __init__(self, cct, components):

        self.has_ic = False
        self.zeroic = True
        self.has_s = False
        self.has_transient = False
        self.ac_count = 0
        self.dc_count = 0
        self.causal = True
        self.reactive = False
        self.independent_sources = []
        self.dependent_sources = []
        self.control_sources = []
        self.reactances = []
        self.ics = []

        for eltname, elt in cct.elements.items():
            if elt.need_control_current:
                self.control_sources.append(elt.args[0])
            if elt.has_ic is not None:
                if elt.has_ic:
                    self.has_ic = True
                    self.ics.append(eltname)
                if not elt.zeroic:
                    self.zeroic = False
            if elt.independent_source:
                self.independent_sources.append(eltname)
                if elt.has_s_transient:
                    self.has_s = True
                if elt.has_transient:
                    self.has_transient = True
                if elt.is_ac:
                    self.ac_count += 1
                if elt.is_dc:
                    self.dc_count += 1
                if not elt.is_causal:
                    self.causal = False
            if elt.dependent_source:
                self.dependent_sources.append(eltname)
            if elt.reactive:
                self.reactive = True
                self.reactances.append(eltname)

        num_sources = len(self.independent_sources)

        self.has_dc = self.dc_count > 0
        self.has_ac = self.ac_count > 0

        self.ac = self.ac_count > 0 and \
            (num_sources == self.ac_count) and not self.has_ic
        self.dc = self.dc_count > 0 and \
            (num_sources == self.dc_count) and not self.has_ic
        self.causal = self.causal and self.zeroic
        self.time_domain = not self.reactive and not self.has_s

        self.ivp = self.has_ic
        self.switching = components.switches != []

        if not self.reactive and self.has_ic:
            raise ValueError('Non-reactive component with initial conditions')

    def __repr__(self):

        return str(vars(self))
