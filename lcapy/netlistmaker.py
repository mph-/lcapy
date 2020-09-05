from .netlisthelper import NetlistHelper

class NetlistMaker(NetlistHelper):

    def __init__(self, net):

        self.net = net
    
    def __call__(self, form='horizontal'):

        if form in ('horizontal', 'default'):
            dir = 'right'
        elif form == 'vertical':
            dir = 'down'
        else:
            raise ValueError('Unknown form ' + form)

        n1 = self._node
        n2 = self._node        
        return self.net._net_make(self, n2, n1, dir)    
