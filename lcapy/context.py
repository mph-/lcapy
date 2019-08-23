class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

        
class Context(object):

    def __init__(self):
        self.symbols = AttrDict()
        self.assumptions = {}
        self.previous = None
        # Noise instance identifier
        self.nid = 0
