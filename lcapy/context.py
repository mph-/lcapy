from sympy.assumptions.assume import global_assumptions

class Context(object):

    def __init__(self):
        self.symbols = {}
        self.assumptions = {}
        self.previous = None
        self.nid = 0

    def new(self):

        new_context = Context()
        new_context.symbols.update(self.symbols)
        new_context.assumptions.update(self.assumptions)
        return new_context

    def switch(self):

        global context

        self.previous = context
        context = self
        global_assumptions.clear()
        global_assumptions.update(self.assumptions)

    def restore(self):

        if self.previous is None:
            return

        self.assumptions.update(global_assumptions)
        global_assumptions.clear()
        global_assumptions.update(self.previous.assumptions)


global_context = Context()
context = global_context

        
