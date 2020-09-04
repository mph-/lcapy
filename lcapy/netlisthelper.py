class NetlistHelper(object):

    @property 
    def _node(self):

        if not hasattr(self, '_node_counter'):
            self._node_counter = 0
        ret = self._node_counter
        self._node_counter += 1
        return ret

    def _make_id(self, kind):
        """Make identifier"""

        if not hasattr(self, '_anon'):
            self._anon = {}            

        if kind not in self._anon:
            self._anon[kind] = 0
        self._anon[kind] += 1
        return self._anon[kind]
    
