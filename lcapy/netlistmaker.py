"""This module provides the NetlistMaker class for generating
horizontal or vertical netlists from a network.

Copyright 2020 Michael Hayes, UCECE

"""

from .netlisthelper import NetlistHelper

class NetlistMaker(NetlistHelper):

    def __init__(self, net, form='horizontal', evalf=None):

        self.net = net
        self.evalf = evalf        

        if form == 'horizontal':
            self.dir = 'right'
        elif form == 'vertical':
            self.dir = 'down'
        else:
            raise ValueError('Unknown form ' + form)

    def __call__(self):

        n1 = self._node
        n2 = self._node        
        return self.net._net_make(self, n2, n1, self.dir)    
