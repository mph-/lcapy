class Assumptions(dict):

    def __init__(self, *args, **kwargs):

        super (Assumptions, self).__init__(*args)

        for assumption, value in kwargs.items():
            self.set(assumption, value)
    
    def set(self, assumption, value):

        if assumption in ('dc', 'ac', 'causal', 'unknown'):
            if value:
                self.pop('dc', None)
                self.pop('ac', None)
                self.pop('causal', None)
                self.pop('unknown', None)                
                self[assumption] = value
            else:
                self.pop(assumption, None)
        else:
            self[assumption] = value

    def get(self, assumption):

        return self[assumption]

    def merge(self, **assumptions):

        for assumption, value in assumptions.items():
            self.set(assumption, value)

    def copy(self):

        return super(Assumptions, self).copy()
    
    def sympy_assumptions(self):
        """Return dict of the SymPy assumptions such as complex, positive, etc."""

        assumptions = {}
        for assumption, value in self.items():
            if assumption not in ('nid', 'ac', 'dc', 'causal', 'unknown'):
                assumptions[assumption] = value
        return assumptions
