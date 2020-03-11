from __future__ import print_function

def unique(alist):

    # Order preserving...  list(set(alist)) gives different results
    # for different runs.
    used = set()
    return [x for x in alist if x not in used and (used.add(x) or True)]


class Cnodes(dict):
    """Common nodes"""

    def __init__(self, nodes):

        super (Cnodes, self).__init__()
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


class Gedge(object):

    def __init__(self, cpt, from_gnode, to_gnode, size, stretch=False):

        self.cpt = cpt
        self.from_gnode = from_gnode
        self.to_gnode = to_gnode
        self.size = size
        self.stretch = stretch
        self.stretched = False

    def __repr__(self):

        return '%s --%s%s--> %s' % (
            self.from_gnode.fmt_name, self.size, '*' if self.stretch else '',
            self.to_gnode.fmt_name)
    
    @property
    def name(self):
        if self.cpt is None:
            return 'dummy'
        return self.cpt.name


class Gnode(object):

    def __init__(self, name):

        self.name = name
        self.prev = None
        self.next = None
        self.dist = 0
        self.pos = None
        self.fedges = []
        self.redges = []

    def add_fedge(self, edge):

        for edge1 in self.fedges:
            if (edge1.cpt == edge.cpt and 
                edge1.to_gnode == edge.to_gnode and
                edge1.from_gnode == edge.from_gnode):
                return

        self.fedges.append(edge)

    def add_redge(self, edge):

        for edge1 in self.redges:
            if (edge1.cpt == edge.cpt and 
                edge1.to_gnode == edge.to_gnode and
                edge1.from_gnode == edge.from_gnode):
                return

        self.redges.append(edge)

    @property
    def fmt_name(self):

        if isinstance(self.name, tuple):
            return ', '.join(self.name)
        return self.name

    def __repr__(self):
        
        return self.fmt_name

    
class Graph(dict):

    def __init__(self, name, nodes, debug=False):

        self.cnodes = Cnodes(nodes)
        self.name = name
        self.debug = debug

    def __repr__(self):

        s = ''
        for gnode in self.values():
            if gnode.fedges == []:
                s += '@%s %s\n' % (gnode.dist, gnode.fmt_name)
            else:
                for edge in gnode.fedges:
                    s += '@%s %s\n' % (gnode.dist, edge)
        return s

    def __getitem__(self, key):

        # If key is an integer, convert to a string.
        if isinstance(key, int):
            key = '%d' % key

        # Allow indexing by a node name rather than by cnode tuple.
        if key in self.cnodes:
            key = self.cnodes[key]

        return super(Graph, self).__getitem__(key)

    def link(self, n1, n2):
        """Make nodes n1 and n2 share common node"""
        
        self.cnodes.link(n1, n2)

    def add(self, cpt, n1, n2, size, stretch):
        """Add cpt between nodes n1 and n2 to the graph"""
        
        if size == 0:
            return

        if size < 0:
            n1, n2 = n2, n1
            size = -size

        if n1 in self.cnodes:
            n1 = self.cnodes[n1]
        if n2 in self.cnodes:
            n2 = self.cnodes[n2]

        gnode1 = self.add_node(n1)
        gnode2 = self.add_node(n2)

        self.add_edges(cpt, gnode1, gnode2, size, stretch)

    def add_node(self, n):
        
        if n not in self:
            self[n] = Gnode(n)
        return self[n]

    def add_edges(self, cpt, gnode1, gnode2, size, stretch):

        gnode1.add_fedge(Gedge(cpt, gnode1, gnode2, size, stretch))
        gnode2.add_redge(Gedge(cpt, gnode2, gnode1, size, stretch))

    @property
    def nodes(self):
        return self.keys()

    @property
    def all_nodes(self):
        return unique(self.cnodes.values())

    def add_start_nodes(self):

        if 'start' in self:
            return

        gnodes = self.values()
        start = self.add_node('start')
        end = self.add_node('end')

        for gnode in gnodes:
            if gnode.redges == []:
                self.add_edges(None, start, gnode, 0, True)
            if gnode.fedges == []:
                self.add_edges(None, gnode, end, 0, True)

    def assign_fixed1(self, gnode):

        for edge in gnode.fedges:
            if (not edge.stretch and edge.to_gnode.pos is not None
                and edge.to_gnode.name != 'end'):
                gnode.pos = edge.to_gnode.pos - edge.size
                return True
        for edge in gnode.redges:
            if (not edge.stretch and edge.to_gnode.pos is not None
                and edge.to_gnode.name != 'start'):
                gnode.pos = edge.to_gnode.pos + edge.size
                return True
        return False

    def assign_fixed(self, unknown):
        """ Assign node positions to nodes with fixed edge lengths to
        nodes with known positions.  Iterate until no more changes.
        This stage is not needed but provides a minor optimisation."""
        
        changes = True
        while changes and unknown != []:
            for n in unknown:
                gnode = self[n]
                changes = self.assign_fixed1(gnode)
                if changes:
                    unknown.remove(n)
                    break

    def assign_stretch1(self, gnode):

        # For each node, find the longest path between nodes with
        # known positions.  The path between the nodes with known
        # positions is then traversed, measuring the total separation
        # and the number of stretchable components.  (Note, there may
        # be multiple paths of the same length; it does not matter
        # which one is chosen.) The stretch is then calculated from
        # the distance of the path, the total separation, and the
        # number of stretches.  If only one known position is found,
        # the node is dangling and the stretch is zero.
        
        # TODO, check for multiple paths with a conflict, for example,
        # a stretchy path of a longser distance than a fixed path.

        to_stretches = 0
        to_dist = 0
        self.longest_path_to_known(gnode)
        to_gnode = gnode
        while to_gnode.next is not None:
            if to_gnode.next.stretch:
                to_stretches += 1
            to_dist += to_gnode.next.size
            to_gnode = to_gnode.next.to_gnode

        from_stretches = 0
        from_dist = 0
        self.longest_path_to_known(gnode, False)
        from_gnode = gnode
        while from_gnode.next is not None:
            if from_gnode.next.stretch:
                from_stretches += 1
            from_dist += from_gnode.next.size
            from_gnode = from_gnode.next.to_gnode

        if from_gnode.name == 'start' and to_gnode.name == 'end':
            # There is a chance that we have naively picked an
            # unlucky gnode. 
            return False

        if from_gnode.name == 'start':
            # Have dangling node, so no stretch needed.
            gnode.pos = to_gnode.pos - to_dist
            return True

        if to_gnode.name == 'end':
            # Have dangling node, so no stretch needed.
            gnode.pos = from_gnode.pos + from_dist
            return True

        separation = to_gnode.pos - from_gnode.pos
        extent = to_dist + from_dist
        if extent - separation > 1e-6:
            raise ValueError('Inconsistent %s schematic graph, component will not fit:  separation %s between %s and %s, need %s.\n%s' % (self.name, separation, from_gnode, to_gnode, extent, self))
            
        if from_stretches == 0:
            gnode.pos = from_gnode.pos + from_dist
        elif to_stretches == 0:
            gnode.pos = to_gnode.pos - to_dist
        else:
            stretch = (separation - extent) / (to_stretches + from_stretches)
            gnode.pos = from_gnode.pos + from_dist + stretch * from_stretches
        return True

    def assign_stretch(self, unknown):
        """Use a worklist algorithm to assign nodes with unknown positions
        that are connected via stretchable edges to the known nodes.

        """

        changes = True
        while changes and unknown != []:
            for n in unknown:
                gnode = self[n]
                changes = self.assign_stretch1(gnode)
                if changes:
                    unknown.remove(n)
                    self.assign_fixed(unknown)
                    break

    def analyse(self, stage=None):

        self.add_start_nodes()

        for gnode in self.values():
            gnode.pos = None

        unknown = list(self.keys())
        unknown.remove('start')
        unknown.remove('end')

        if unknown == []:
            pos = {}
            for n, gnode in self.cnodes.items():
                pos[n] = 0
            return pos, 0

        # Find longest path through the graph.
        self.longest_path(self['start'])

        try:
            # Nodes on the longest path have known positions.
            gnode = self['end']
            while gnode != None and gnode.name != 'start':
                if gnode.name in unknown:
                    unknown.remove(gnode.name)
                gnode.path = True
                gnode.pos = gnode.dist
                gnode.fixed = True
                gnode = gnode.prev.from_gnode
                self['start'].pos = 0
        except AttributeError:
            msg = "The %s schematic graph is dodgy, probably a "
            "component is attached to the wrong nodes.\n" % self.name
            if self.debug:
                msg += str(self)
            raise AttributeError(msg)

        if stage == 1:
            return
            
        self.assign_fixed(unknown)

        if stage == 2:
            return

        self.assign_stretch(unknown)

        if unknown != []:
            raise ValueError('Cannot assign nodes %s for %s graph.\n%s' %
                             (unknown, self.name, self))
            
        self.check_positions()
  
        try:
            pos = {}
            for n, gnode in self.cnodes.items():
                pos[n] = self[gnode].pos

        except KeyError:
            # TODO determine which components are not connected.
            msg = "The %s schematic graph is dodgy, probably a component is attached to the wrong nodes.\n" % self.name
            if self.debug:
                msg += str(self)
            raise AttributeError(msg)            

        distance_max = self['end'].pos

        return pos, distance_max

    def longest_path_to_known(self, start, forward=True):
        """Find longest path through DAG to a node with a known dist."""

        def traverse(gnode):

            if gnode.name in ('start', 'end'):
                # Choose as last resort
                gnode.next = None
                return 1000

            if gnode.pos is not None:
                gnode.next = None
                return gnode.pos

            edges = gnode.fedges if forward else gnode.redges
            min_dist = 2000
            for edge in edges:
                next_gnode = edge.to_gnode
                # TODO, memoize
                dist = traverse(next_gnode) - edge.size
                if dist < min_dist:
                    min_dist = dist
                    next_gnode.prev = edge
                    gnode.next = edge
            return min_dist

        start.dist = 0
        try:
            return traverse(start)
        except RuntimeError:
            msg = "The %s schematic graph is dodgy, probably a component is attached to the wrong nodes.\n" % self.name
            if self.debug:
                msg += str(self)
            raise RuntimeError(msg)


    def longest_path(self, start, forward=True):
        """Find longest path through DAG."""

        for gnode in self.values():
            gnode.dist = -1
            gnode.prev = None
            gnode.next = None

        def traverse(gnode):

            if gnode.pos is not None:
                return True

            edges = gnode.fedges if forward else gnode.redges
            for edge in edges:
                next_gnode = edge.to_gnode
                dist = gnode.dist + edge.size
                if dist > next_gnode.dist:
                    next_gnode.dist = dist
                    next_gnode.prev = edge
                    gnode.next = edge
                    if traverse(next_gnode):
                        return True
            return False

        start.dist = 0
        try:
            traverse(start)
        except RuntimeError:
            msg = "The %s schematic graph is dodgy, probably a component is attached to the wrong nodes.\n" % self.name
            if self.debug:
                msg += str(self)
            
            raise RuntimeError(msg)

    def check_positions(self):

        for gnode in self.values():
            for edge in gnode.fedges:
                dist = edge.to_gnode.pos - gnode.pos
                if edge.stretch:
                    if dist - edge.size < -1e-6:
                        print('Distance conflict %s vs %s in %s graph for %s between nodes %s and %s,'
                              ' due to incompatible sizes' % (
                                  dist, edge.size, self.name,
                                  edge.name, gnode.fmt_name,
                                  edge.to_gnode.fmt_name))
                else:
                    if abs(dist - edge.size) > 1e-6:
                        print('Stretch conflict %s vs %s in %s graph for %s between nodes %s and %s,'
                              ' due to incompatible sizes' % (
                                  dist, edge.size, self.name,
                                  edge.name, gnode.fmt_name,
                                  edge.to_gnode.fmt_name))

    def dot(self, filename=None, stage=None):
        """Generate directed graph using graphviz notation"""

        from .system import tmpfilename, run_dot
        from .schematic import display_matplotlib
        from os import path
        
        def fmt_dec(value):
            return ('%.2f' % value).rstrip('0').rstrip('.')

        if filename is None:
            filename = tmpfilename('.png')
            self.dot(filename=filename, stage=stage)
            display_matplotlib(filename)
            return

        base, ext = path.splitext(filename)
        if ext in ('.pdf', '.png'):
            dotfilename = filename + '.dot'
            self.dot(dotfilename, stage=stage)
            run_dot(dotfilename, filename)
            return

        if stage != 0:
            self.analyse(stage=stage)

        dotfile = open(filename, 'w')
        dotfile.write('strict digraph {\n\tgraph [rankdir=LR];\n')

        for gnode in self.values():
            if hasattr(gnode, 'fixed'):
                colour = 'yellow'
            else:
                colour = 'red' if gnode.pos is not None else 'blue'
            if gnode.name in ('start', 'end'):
                colour = 'green'

            pos = gnode.pos
            if pos is None or pos < 1e-6:
                pos = 0

            dotfile.write('\t"%s"\t [style=filled, color=%s, xlabel="@%s"];\n' % (gnode.fmt_name, colour, fmt_dec(pos)))

        for gnode in self.values():
            for edge in gnode.fedges:
                colour = 'black' if edge.stretch else 'red'
                dotfile.write('\t"%s" ->\t"%s" [ color=%s, label="%s%s" ];\n' % (
                    gnode.fmt_name, edge.to_gnode.fmt_name, colour, 
                    fmt_dec(edge.size), '*' if edge.stretch else ''))

        dotfile.write('}\n')
        dotfile.close()

