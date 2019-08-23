from sympy.assumptions.assume import global_assumptions

from .context import Context

class State(object):

    def __init__(self):
        self.global_context = Context()
        self.context = self.global_context
        self.previous_context = None

    def switch_context(self, context):
        
        self.previous_context = self.context
        self.context = context

        global_assumptions.clear()
        global_assumptions.update(self.context.assumptions)

        
    def restore_context(self):

        self.context.assumptions.update(global_assumptions)
        
        self.context = self.previous_context
        self.previous_context = None
        
        global_assumptions.clear()
        global_assumptions.update(self.context.assumptions)

state = State()
