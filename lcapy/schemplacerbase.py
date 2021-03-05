from .schemmisc import Pos


class SchemPlacerBase(object):



    def solve(self, node_spacing):

        self._make_graphs()

        xpos, width = self.xgraph.solve()
        ypos, height = self.ygraph.solve()

        scale = node_spacing
        for n, node in self.nodes.items():
            node.pos = Pos(xpos[n] * scale, ypos[n] * scale)

        return width, height
        
