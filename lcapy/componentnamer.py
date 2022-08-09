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

    def name(self, cpt_type, names=None):

        cpt_id = self.cpt_id(cpt_type)

        newname = cpt_type + cpt_id
        # Ensure that name is unique if list of names supplied.
        if names is not None and newname in names:
            return self.name(cpt_type, names)

        return newname

    def __call__(self, cpt_type, names=None):
        """Create component name based on `cpt_type` (`R`, `L`, etc.)
        and the name_kind (`anon`, `t`, etc).

        `names` is a list of known names to avoid."""

        return self.name(cpt_type, names)
