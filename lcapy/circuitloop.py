class CircuitEdge(list):

    def __init__(self, node1, node2, cpt):

        self.node1 = node1
        self.node2 = node2
        self.cpt = cpt

    @property
    def name(self):

        if self.cpt is None:
            return ''
        return self.cpt.name

    def __repr__(self):

        return f'CircuitEdge({self.node1}, {self.node2}, {self.name})'

    def has_cpt_name(self, cpt_name):

        if self.cpt is None:
            return False
        return self.cpt.name == cpt_name


class CircuitLoop(list):

    def __repr__(self):

        return 'CircuitLoop((' + ', '.join([repr(edge) for edge in
                                            self]) + '))'

    @classmethod
    def from_nodes_cpts(cls, nodes, cpts):

        edges = []

        if len(nodes) != len(cpts):
            raise ValueError('Wrong number of cpts for nodes')

        for m in range(len(cpts) - 1):
            edges.append(CircuitEdge(nodes[m], nodes[m + 1], cpts[m]))
        edges.append(CircuitEdge(nodes[-1], nodes[0], cpts[-1]))

        return cls(edges)

    def has_cpt(self, cpt):

        for edge in self:
            if edge.cpt is None:
                continue
            if edge.cpt == cpt:
                return True
        return False

    def has_cpt_name(self, cpt_name):

        for edge in self:
            if edge.cpt is None:
                continue
            if edge.cpt.name == cpt_name:
                return True
        return False

    def edge_by_cpt_name(self, cpt_name):

        for edge in self:
            if edge.cpt is None:
                continue
            if edge.cpt.name == cpt_name:
                return edge
        raise ValueError('Unknown cpt ' + cpt_name)

    def as_nodes(self):

        nodes = []
        for edge in self:
            nodes.append(edge.node1)
        return nodes

    def as_cpts(self):

        cpts = []
        for edge in self:
            if edge.cpt is not None:
                cpts.append(edge.cpt)
        return cpts

    def as_cpt_names(self):

        cpt_names = []
        for edge in self:
            if edge.cpt is not None:
                cpt_names.append(edge.cpt.name)
        return cpt_names
