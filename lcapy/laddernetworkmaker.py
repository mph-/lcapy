from lcapy import Circuit
from lcapy.circuitgraph import Edges, CircuitGraph


class LadderNetworkMaker:
    """This class converts a circuit into an unbalanced ladder network."""

    def __init__(self, cct, debug=False):

        self.cct = cct
        self.cg = CircuitGraph.from_circuit(cct)
        self.debug = debug

    def _find(self, N1p, N1m, N2p, N2m):
        """Convert a circuit into an unbalanced ladder network.
        The input port is defined by the nodes `N1p` and `N1m`.
        The output port is defined by the nodes `N2p` and `N2m`.

        The nodes `N1p` and `N1m` must be the same.

        This returns `None` if the network does not have a ladder
        topology.  Otherwise, it returns a list of alternating series
        and shunt paths, where a path is one or more edges.  If there
        is no initial series path, it is assigned `[]`.

        """

        # Rename nodes to the unique nodes used in the circuitgraph
        N1p, N1m, N2p, N2m = self.cg.canonical_nodes(N1p, N1m, N2p, N2m)

        parts = []

        if N1m != N2m:
            return parts

        cg = self.cg
        node = N1p
        initial = True

        while True:

            # Can have multiple series paths to the next node
            # Let's assume only a single path
            opts = []

            edges = cg.node_edges(node)

            if len(edges) != 1:
                if initial:
                    # No initial series branch
                    parts.append([])
                else:
                    return []
            else:
                edge = edges[0]

                path = self.cg.series_path(edge, None, N2p)
                if path == []:
                    return []
                if self.debug:
                    print('series: ' + str(path))
                opts.append(path)
                self.cg.remove_edges(path)
                edge = path[-1]
                node = edge.to_node
                parts.append(opts)

            initial = False

            # _Find all the parallel paths
            opts = []
            edges = cg.node_edges(node)
            for edge in edges:
                path = self.cg.series_path(edge, N2m, N2p)
                if path != []:
                    if self.debug:
                        print('parallel: ' + str(path))
                    opts.append(path)
                    self.cg.remove_edges(path)
            parts.append(opts)

            if node == N2p:
                return parts

    def make(self, N1p, N1m, N2p, N2m):
        """Return two-port unbalanced ladder network or `None` if the netlist
        does not have a ladder topology between the specified nodes.

        The input port is defined by the nodes `N1p` and `N1m`.
        The output port is defined by the nodes `N2p` and `N2m`.

        The nodes `N1p` and `N1m` must be the same."""

        if self.cct.dependent_sources != []:
            return None
        if self.cct.mutual_couplings != []:
            return None
        if self.cct.transformers != []:
            return None

        foo = self._find(N1p, N1m, N2p, N2m)

        if foo == []:
            return None

        if foo[-1] == []:
            foo = foo[0:-1]

        from lcapy.oneport import Ser, Par
        from lcapy.twoport import Ladder, LadderAlt

        args = []
        for p in foo:
            pars = []
            for q in p:
                cpts = [
                    self.cct[q1.cpt_name].cpt for q1 in q if q1.cpt_name[0] != 'W']
                if len(cpts) != 1:
                    cpt = Ser(*cpts)
                else:
                    cpt = cpts[0]
                pars.append(cpt)

            if len(pars) == 0:
                cpt = None
            elif len(pars) == 1:
                cpt = pars[0]
            else:
                cpt = Par(*pars)

            args.append(cpt)

        if args[0] is None:
            return LadderAlt(*args[1:])

        return Ladder(*args)


def test():

    cct = Circuit("""
P1  1   0   ; down
R_1 1   2   ; right
C_1 2   3   ; right
C_2 3   4   ; down
L_1 4   0_2 ; down
C_3 3   5   ; right
C_4 5   6   ; down
L_2 6   0_3 ; down
C_5 5   7   ; right
C_7 7   8   ; down
L_3 8   0_4 ; down
C_6 7   9   ; right
R_2 9   0_5 ; down
P2  9_a 0_6 ; down
W   0   0_2 ; right
W   0_2 0_3 ; right
W   0_3 0_4 ; right
W   0_4 0_5 ; right
W   0_5 0_6 ; right
W   9   9_a ; right
""")

    lm = LadderNetworkMaker(cct)
    tp = lm.make(1, 0, 9, 0)
    print(tp)
    return tp
