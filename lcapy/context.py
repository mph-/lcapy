"""This module provides the Context class that maintains the context
for a Circuit.

Copyright 2014--2020 Michael Hayes, UCECE

"""

from .attrdict import AttrDict

class Context(object):

    def __init__(self):
        self.symbols = AttrDict()
        self.assumptions = {}
        # Noise instance identifier
        self.nid = 0
