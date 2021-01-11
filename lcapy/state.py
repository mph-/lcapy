"""This module provides the State class that maintains the global state,
such as the symbol context.

Copyright 2014--2021 Michael Hayes, UCECE

"""

from sympy.assumptions.assume import global_assumptions

from .context import Context
from .config import use_units
from copy import copy

class State(object):

    def __init__(self):
        self.global_context = Context()
        self.context = self.global_context
        self.previous_context = []
        self.use_units = use_units

    def new_context(self):

        context = Context()
        context.symbols = copy(self.global_context.symbols)
        return context
        
    def switch_context(self, context):

        self.previous_context.append(self.context)
        self.context = context

        global_assumptions.clear()
        global_assumptions.update(self.context.assumptions)

    def restore_context(self):

        self.context.assumptions.update(global_assumptions)
        
        self.context = self.previous_context.pop()
        
        global_assumptions.clear()
        global_assumptions.update(self.context.assumptions)

state = State()
