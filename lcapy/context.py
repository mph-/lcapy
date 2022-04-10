"""This module provides the Context class that maintains the context
for a Circuit.

Copyright 2014--2022 Michael Hayes, UCECE

"""


class Context(object):

    def __init__(self):
        self.assumptions = {}
        # Noise instance identifier
        self.nid = 0
