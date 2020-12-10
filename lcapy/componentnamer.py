"""This module provides the ComponentNamer class for automatic
naming of components.

Copyright 2020 Michael Hayes, UCECE

"""


class ComponentNamer(object):

    def __init__(self):

        self.names = {}

    def make(self, cpt_kind, name_kind='anon'):

        if name_kind not in self.names:
            self.names[name_kind] = {}

        foo = self.names[name_kind]
            
        if cpt_kind not in foo:
            foo[cpt_kind] = 0

        foo[cpt_kind] += 1

        return cpt_kind + name_kind + str(foo[cpt_kind])
    
    def __call__(self, cpt_kind, name_kind='anon'):
        """Create component name based on `cpt_kind` (`R`, `L`, etc.)
        and the name_kind (`anon`, `t`, etc)."""

        return self.make(cpt_kind, name_kind)
