"""This module provides classes for schematic layout.

Copyright 2014--2020 Michael Hayes, UCECE

"""


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
    """Edge between common nodes"""

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
    """Drawing node"""

    def __init__(self, name):

        self.name = name
        self.dist = 0
        self._pos = None
        self.fedges = []
        self.redges = []

    @property
    def pos(self):
        return self._pos

    @pos.setter    
    def pos(self, value):

        if value is not None and self._pos is not None:
            raise RuntimeError('Changing node %s pos' % self)
        if False and value is not None:
            print('Setting node %s to pos %f' % (self, value))
            
        self._pos = value

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
            return '(' + ', '.join(self.name) + ')'
        return self.name

    def __repr__(self):
        
        return self.fmt_name

    
class GraphPath(list):
    """This is a list of edges defining a path through a graph."""

    @property
    def dist(self):

        distance = 0
        for edge in self:
            distance += edge.size
        return distance

    @property
    def stretches(self):

        stretches = 0
        for edge in self:
            if edge.stretch:
                stretches += 1
        return stretches

    def on(self, match_node):

        if self == []:
            return False
        
        for edge in self:
            if edge.from_gnode == match_node:
                return True

        return edge.to_gnode == match_node

    def to_dist(self, match_node):

        if match_node == self.from_gnode:
            return 0
        
        distance = 0
        for edge in self:
            distance += edge.size
            if edge.to_gnode == match_node:
                return distance
        raise ValueError('Node not on path')

    def from_dist(self, match_node):

        if match_node == self.to_gnode:
            return 0        

        distance = 0
        for edge in reversed(self):
            distance += edge.size
            if edge.from_gnode == match_node:
                return distance
        raise ValueError('Node not on path')    

    @property    
    def gnodes(self):
        """Return list of nodes along path including end nodes."""
        
        nodes = []
        nodes.append(self.from_gnode)
        for edge in self:
            nodes.append(edge.to_gnode)
        return nodes
    
    @property
    def to_gnode(self):

        return self[-1].to_gnode

    @property
    def from_gnode(self):

        return self[0].from_gnode

                
class Graph(dict):
    """Drawing graph"""

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
        """Nodes without forward edges are connected to the end node.
        Nodes without reverse edges are connected to the start node."""

        if 'start' in self:
            return

        start = self.add_node('start')
        end = self.add_node('end')

        gnodes = self.values()
        for gnode in gnodes:
            if gnode.redges == [] and gnode != start:
                self.add_edges(None, start, gnode, 0, True)
            if gnode.fedges == [] and gnode != end:
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

    def assign_stretchy1(self, gnode, unknown):

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
        # a stretchy path of a longer distance than a fixed path.
        
        to_path = self.path_to_closest_known(gnode, forward=True)
        from_path = self.path_to_closest_known(gnode, forward=False)

        to_gnode = to_path.to_gnode
        # Note, from_path is reversed (the edges are backward edges)
        from_gnode = from_path.to_gnode

        if from_gnode.name == 'start' and to_gnode.name == 'end':
            # There is a chance that we have naively picked an
            # unlucky gnode. 
            return False

        if from_gnode.name == 'start' and len(to_path) == 1:
            # Have dangling node, so no stretch needed.
            gnode.pos = to_gnode.pos - to_path.dist
            unknown.remove(gnode.name)
            return True

        if to_gnode.name == 'end' and len(from_path) == 1:
            # Have dangling node, so no stretch needed.
            gnode.pos = from_gnode.pos + from_path.dist
            unknown.remove(gnode.name)
            return True

        path = self.longest_path(from_gnode, to_gnode)

        if len(path) == 1:
            return False
        
        # The gnode may not be on this path but this does not matter.
        # It will be processed asain later.

        stretches = path.stretches
        separation = to_gnode.pos - from_gnode.pos
        extent = path.dist        

        if extent - separation > 1e-6:
            raise ValueError('Inconsistent %s schematic graph, component(s) will not fit:  separation %s between %s and %s, need %s.\n%s' % (self.name, separation, from_gnode, to_gnode, extent, self))

        
        # This how much each component needs to stretch.
        stretch = (separation - extent) / stretches

        pos = from_gnode.pos
        for edge in path:
            pos += edge.size
            if edge.stretch:
                pos += stretch

            if edge != path[-1]:
                edge.to_gnode.pos = pos
                unknown.remove(edge.to_gnode.name)

        return True

    def assign_stretchy(self, unknown):
        """Use a worklist algorithm to assign nodes with unknown positions
        that are connected via stretchable edges to the known nodes."""

        changes = True
        while changes and unknown != []:
            for n in unknown:
                gnode = self[n]
                if gnode.pos is not None:
                    unknown.remove(n)
                    break
                
                changes = self.assign_stretchy1(gnode, unknown)
                if changes:
                    self.assign_fixed(unknown)
                    break

    def assign_longest(self, path, unknown):

        if path == []:
            raise RuntimeError('Empty path')
        
        pos = 0
        for edge in path:
            edge.from_gnode.pos = pos
            pos += edge.size
            unknown.remove(edge.from_gnode.name)

        edge = path[-1]
        edge.to_gnode.pos = pos
        unknown.remove(edge.to_gnode.name)

    def prune(self):
        """Remove redundant paths from graph.  If there are two
        parallel stretchy edges, the longest is chosen."""

        for gnode in self.values():

            def check(edges, grizzle=False):

                # If don't have multiple edges, cannot have redundancy
                if len(edges) < 2:
                    return edges
                
                fromto = {}
                for edge in edges:
                    key = (edge.from_gnode, edge.to_gnode)
                    if key not in fromto:
                        fromto[key] = []
                    fromto[key].append(edge)

                newedges = []                        
                for key, edges in fromto.items():

                    best_edge = edges[0]
                    fixed = False
                    size = 0
                    for edge in edges:
                        if not edge.stretch:
                            fixed = True
                            size = edge.size
                            best_edge = edge
                            break

                    if fixed and grizzle:
                        for edge in edges:
                            if edge.size != size and not edge.stretch:
                                print('Component fixed size violation for %s, expected %f, got %f' % (edge.cpt, size, edge.size))
                            elif edge.size > size:                                
                                print('Component size violation for %s, expected %f or less, got %f' % (edge.cpt, size, edge.size))
                                
                    elif not fixed:
                        # Choose largest component size since all are stretchy
                        for edge in edges:
                            if edge.size > size:
                                size = edge.size
                                best_edge = edge

                    newedges.append(best_edge)

                return newedges
                    

            gnode.fedges = check(gnode.fedges, grizzle=True)
            gnode.redges = check(gnode.redges)            

    def analyse(self, stage=None):
        """
        Analyse graph assigning gnode positions.

        stage 0 --- add start/end gnodes
        stage 1 --- prune redundant edges
        stage 2 --- assign gnode positions on longest path
        stage 3 --- assign gnode positions with fixed positions to known gnodes
        stage 4 --- assign all gnode positions
        """
        
        unknown = list(self.keys())

        if unknown == []:
            # If there are no constraints, all the nodes have distance 0.
            pos = {}
            for n, gnode in self.cnodes.items():
                pos[n] = 0
            return pos, 0
        
        self.add_start_nodes()
        unknown.append('start')
        unknown.append('end')        

        if stage == 0:
            return

        # Prune redundant edges from graph.        
        self.prune()

        if stage == 1:
            return
            
        # Find longest path through the graph.  This provides the
        # dimension for the graph.
        path = self.longest_path(self['start'], self['end'])

        # Assign gnodes on the longest path.        
        self.assign_longest(path, unknown)
        
        if stage == 2:
            return

        # Assign gnodes at a fixed distance from nodes with known positions.
        self.assign_fixed(unknown)

        if stage == 3:
            return

        # Assign remaining gnodes.
        self.assign_stretchy(unknown)

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

    def longest_path(self, from_gnode, to_gnode, forward=True):
        """Find longest path through DAG from from_gnode to to_gnode
        or to a gnode with a known position."""

        # gnode.dist specifies the distance to the to_gnode.  It is
        # None if the distance has yet to be determined and is -1 if
        # the to_gnode cannot be reached from gnode.
        for gnode in self.values():
            gnode.dist = None
            gnode.next = None            

        def traverse(gnode):

            # If have visited this gnode before then know dist.
            if gnode.dist is not None:
                return gnode.dist

            # If have reached goal, return dist of 0.
            if gnode == to_gnode:
                gnode.dist = 0                
                return 0

            # Do not traverse nodes at known positions.
            if gnode.pos is not None and gnode != from_gnode:
                return -1

            gnode.dist = -1

            edges = gnode.fedges if forward else gnode.redges
            for edge in edges:
                dist = traverse(edge.to_gnode)
                if dist >= 0 and gnode.dist < dist + edge.size:
                    gnode.dist = dist + edge.size
                    # Mark current best edge
                    gnode.next = edge

            return gnode.dist

        try:
            traverse(from_gnode)
        except RuntimeError:
            msg = "The %s schematic graph is dodgy, probably a component is attached to the wrong nodes.\n" % self.name
            if self.debug:
                msg += str(self)
            
            raise RuntimeError(msg)

        return self.makepath(from_gnode)

    def makepath(self, from_gnode):
        
        # Construct path
        path = GraphPath()
        gnode = from_gnode
        while gnode.next is not None:
            path.append(gnode.next)
            gnode = gnode.next.to_gnode
            
        return path
        
    def path_to_closest_known(self, from_gnode, forward=True):
        """Find path through DAG to the closest node with a known pos."""

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

        from_gnode.dist = 0
        try:
            traverse(from_gnode)
        except RuntimeError:
            # W a b; right
            # W c d; right
            # W a c; down
            # W b d; up
            msg = "The %s schematic graph has a loop.  For example, a node needs to be both above and below another node.  Probably a component is attached to the wrong nodes.\n" % self.name
            if self.debug:
                msg += str(self)
            raise RuntimeError(msg)

        return self.makepath(from_gnode)
        
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

        self.analyse(stage=stage)

        dotfile = open(filename, 'w')
        dotfile.write('digraph {\n\tgraph [rankdir=LR];\n')

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

