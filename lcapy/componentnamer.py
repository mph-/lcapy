"""This module provides the ComponentNamer class for automatic
naming of components.

Copyright 2020--2022 Michael Hayes, UCECE

"""


class ComponentNamer(object):

    def __init__(self):

        self.counts = {}

    def cpt_id(self, cpt_type):

        if cpt_type not in self.counts:
            self.counts[cpt_type] = 0

        self.counts[cpt_type] += 1

        return str(self.counts[cpt_type])

    def name(self, cpt_type):

        cpt_id = self.cpt_id(cpt_type)

        return cpt_type + cpt_id

    def __call__(self, cpt_type):
        """Create component name based on `cpt_type` (`R`, `L`, etc.)
        and the name_kind (`anon`, `t`, etc)."""

        return self.name(cpt_type)
