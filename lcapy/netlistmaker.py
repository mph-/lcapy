from .netlisthelper import NetlistHelper

class NetlistMaker(NetlistHelper):

    def __init__(self, net):

        self.net = net
    
    def __call__(self):
    
        return self.net._net_make(self)    
