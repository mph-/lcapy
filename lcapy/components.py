"""
This module collects collates the circuit components.

Copyright 2022 Michael Hayes, UCECE

"""


class Components:

    def __init__(self, cct):

        self.capacitors = []
        self.current_sources = []
        self.inductors = []
        self.resistors = []
        self.switches = []
        self.open_circuits = []
        self.ports = []
        self.transformers = []
        self.voltage_sources = []
        self.wires = []

        for eltname, elt in cct.elements.items():

            if elt.is_capacitor:
                self.capacitors.append(eltname)
            elif elt.is_current_source:
                self.current_sources.append(eltname)
            elif elt.is_inductor:
                self.inductors.append(eltname)
            elif elt.is_open_circuit:
                self.open_circuits.append(eltname)
            elif elt.is_port:
                self.ports.append(eltname)
            elif elt.is_resistor:
                self.resistors.append(eltname)
            elif elt.is_switch:
                self.switches.append(eltname)
            elif elt.is_transformer:
                self.transformers.append(eltname)
            elif elt.is_voltage_source:
                self.voltage_sources.append(eltname)
            elif elt.is_wire:
                self.wires.append(eltname)

    def __repr__(self):

        return str(vars(self))
