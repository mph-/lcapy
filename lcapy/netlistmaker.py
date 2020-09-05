from .netlisthelper import NetlistHelper

class NetlistMaker(NetlistHelper):

    def __init__(self, net):

        self.net = net
    
    def __call__(self):

        n1 = self._node
        n2 = self._node        
        return self.net._net_make(self, n2, n1)    
