"""This module provides the ComponentNamer class for automatic
naming of components.

Copyright 2020--2023 Michael Hayes, UCECE

"""


class ComponentNamer(object):
    """
    Is used to provide a unique naming of components until it is reset.
    """
    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super(ComponentNamer, cls).__new__(cls)
        return cls.instance

    def __init__(self):
        if not hasattr(self, 'names'):
            self.names = []

    def reset(self):
        """Because this is implemented as a singleton a reset can be necessary"""
        self.names = []

    def name(self, cpt_type, names):

        m = 1
        while True:
            name = cpt_type + str(m)
            if name not in names and name not in self.names:
                self.names.append(name)
                return name
            m += 1

    def __call__(self, cpt_type, names=None):
        """Create component name based on `cpt_type` (`R`, `L`, etc.)
        and the name_kind (`anon`, `t`, etc).

        `names` is a list of known names to avoid."""

        return self.name(cpt_type, names)
