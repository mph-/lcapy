def unique(alist):

    # Order preserving...  list(set(alist)) gives different results
    # for different runs.
    used = set()
    return [x for x in alist if x not in used and (used.add(x) or True)]


class Cnodes(dict):
    """Common nodes"""

    def __init__(self, nodes):

        super(Cnodes, self).__init__()
        for node in nodes:
            # Use tuple so hashable.
            self[node] = (node, )

    def link(self, n1, n2):
        """Make nodes n1 and n2 share common node"""

        set1 = self[n1]
        set2 = self[n2]

        newset = tuple(unique(set1 + set2))

        for n in self[n1]:
            self[n] = newset
        for n in self[n2]:
            self[n] = newset

    def all_nodes(self):

        return unique(self.values())
