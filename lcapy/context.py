"""This module provides the Context class that maintains the context
for a Circuit.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from .attrdict import AttrDict

class Context(object):

    def __init__(self):
        self.symbols = AttrDict()
        # User symbols are created with symbol or symbols.
        # These are printed as defined.
        self.user_symbols = {}
        # Domain symbols are symbols such as tsym, ssym.
        self.domain_symbols = {}
        self.assumptions = {}
        # Noise instance identifier
        self.nid = 0
