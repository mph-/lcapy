class Context(object):

    def __init__(self):
        self.symbols = {}
        self.assumptions = {}
        self.previous = None
        # Noise instance identifier
        self.nid = 0
