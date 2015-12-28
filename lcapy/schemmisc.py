import numpy as np

class Pos(object):

    def __init__(self, x, y):

        self.x = x
        self.y = y

    def __mul__(self, scale):

        return Pos(self.x * scale, self.y * scale)

    def __str__(self):

        xstr = ('%.2f' % self.x).rstrip('0').rstrip('.')
        ystr = ('%.2f' % self.y).rstrip('0').rstrip('.')

        return "%s,%s" % (xstr, ystr)

    def __repr__(self):

        return 'Pos(%s)' % self

    @property
    def xy(self):

        return np.array((self.x, self.y))


class Opts(dict):

    def _parse(self, string):

        for part in string.split(','):
            part = part.strip()
            if part == '':
                continue

            if part in ('up', 'down', 'left', 'right'):
                self['dir'] = part
                continue

            fields = part.split('=')
            key = fields[0].strip()
            arg = '='.join(fields[1:]).strip() if len(fields) > 1 else ''
            self[key] = arg

    def __init__(self, arg=None):

        if arg is None:
            return

        if isinstance(arg, str):
            self._parse(arg)
            return

        for key, val in arg.iteritems():
            self[key] = val

    @property
    def size(self):

        size = self.get('size', 1)
        return float(size)

    def format(self):

        return ', '.join(['%s=%s' % (key, val)
                          for key, val in self.iteritems()])

    def copy(self):
        
        return self.__class__(super(Opts, self).copy())

    def strip(self, *args):

        stripped = Opts()
        for opt in args:
            if opt in self:
                stripped[opt] = self.pop(opt)        
        return stripped

    def strip_voltage_labels(self):

        return self.strip('v', 'vr', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<')

    def strip_current_labels(self):

        return self.strip('i', 'ir', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<',
                          'i>_', 'i<_', 'i>^', 'i<^')

    def strip_labels(self):

        return self.strip('l', 'l^', 'l_')

    def strip_all_labels(self):

        self.strip_voltage_labels()
        self.strip_current_labels()
        self.strip_labels()


class EngFormat(object):

    def __init__(self, value, unit=''):

        self.value = value
        self.unit = unit

    def latex(self, trim=True, hundreds=False):

        prefixes = ('f', 'p', 'n', '$\mu$', 'm', '', 'k', 'M', 'G', 'T')

        sfmax = 3

        value = self.value
        m = math.log10(abs(value))

        if m < -1 or m >= 3.0:
            if hundreds:
                # Generate 100 m
                n = int(math.floor(m / 3))
                k = int(math.floor(m)) - n * 3
            else:
                # Generate 0.1
                n = int(round(m / 3))
                k = int(round(m)) - n * 3
        else:
            n = 0
            k = m - 1

        dp = sfmax - k

        idx = n + 5
        if idx < 0:
            idx = 0
            return '%e\,' % value + self.unit
        elif idx >= len(prefixes):
            idx = len(prefixes) - 1
            return '%e\,' % value + self.unit

        fmt = '%%.%df' % dp

        n = idx - 5
        value = value * 10**(-3 * n)

        string = fmt % value

        if trim:
            # Remove trailing zeroes after decimal point
            string = string.rstrip('0').rstrip('.')

        return string + '\,' + r'\mbox{' + prefixes[idx] + self.unit + r'}'


class Cnodes(object):

    """Common nodes"""

    def __init__(self, nodes):

        self.sets = {}
        for node in nodes:
            self.sets[node] = (node, )

    def link(self, n1, n2):
        """Make nodes n1 and n2 share common node"""

        set1 = self.sets[n1]
        set2 = self.sets[n2]
        newset = set1 + set2

        for n in self.sets[n1]:
            self.sets[n] = newset
        for n in self.sets[n2]:
            self.sets[n] = newset

    def _analyse(self):

        # Add dummy cnode at start
        unique = ['dummy'] + list(set(self.sets.values()))
        node_map = {}
        for node, nodes in self.sets.iteritems():
            node_map[node] = unique.index(nodes)

        self._node_map = node_map
        self._nodes = unique

    @property
    def node_map(self):
        """Return mapping of node number to common node number"""

        if not hasattr(self, '_node_map'):
            self._analyse()

        return self._node_map

    def map(self, nodes):

        if not isinstance(nodes, (tuple, list)):
            nodes = list([nodes])

        return [self.node_map[node] for node in nodes]

    @property
    def nodes(self):
        """Return mapping of common node number to tuple of shared nodes"""

        if not hasattr(self, '_nodes'):
            self._analyse()

        return self._nodes

    @property
    def size(self):
        """Return number of common nodes"""

        return len(self.nodes)


class Graph(dict):

    def __init__(self, size, name):

        self.name = name
        for m in range(size):
            self[m] = []

    def add(self, n1, n2, size):
        self[n1].append((n2, size))


    def longest_path(self):
        """Find longest path through DAG.  all_nodes is an iterable for all
        the nodes in the graph, from_nodes is a directory indexed by node
        that stores a tuple of tuples.  The first tuple element is the
        parent node and the second element is the minimium size of the
        component connecting the nodes.
        """

        if self[0] == []:
            raise ValueError("Cannot find start node for graph '%s'. "
                             "Probably a component has an incorrect direction."
                             % self.name)

        all_nodes = self.keys()
        from_nodes = self

        memo = {}

        def get_longest(to_node):

            if to_node in memo:
                return memo[to_node]

            best = 0
            for from_node, size in from_nodes[to_node]:
                best = max(best, get_longest(from_node) + size)

            memo[to_node] = best

            return best

        try:
            length, node = max([(get_longest(to_node), to_node)
                                for to_node in all_nodes])
        except RuntimeError:
            raise RuntimeError(
                ("The schematic graph '%s' is dodgy, probably a component"
                 " is connected to the wrong node\n%s") 
                % (self.name, from_nodes))

        return length, node, memo


class Graphs(object):

    def __init__(self, size, name):

        self.fwd = Graph(size, 'forward ' + name)
        self.rev = Graph(size, 'reverse ' + name)
        self.size = size

    def add(self, n1, n2, size):
        self.fwd.add(n1, n2, size)
        self.rev.add(n2, n1, size)

    @property
    def nodes(self):
        return self.fwd.keys()

    def add_start_nodes(self):

        # Chain all potential start nodes to node 0.
        orphans = []
        rorphans = []
        for m in range(1, self.size):
            if self.fwd[m] == []:
                orphans.append((m, 0))
            if self.rev[m] == []:
                rorphans.append((m, 0))
        self.fwd[0] = rorphans
        self.rev[0] = orphans


class Node(object):

    def __init__(self, name):

        self.name = name
        self._port = False
        self._count = 0
        parts = name.split('_')
        self.rootname = parts[0] if name[0] != '_' else name
        self.primary = len(parts) == 1
        self.list = []

    def append(self, elt):
        """Add new element to the node"""

        if isinstance(elt, P):
            self._port = True

        self.list.append(elt)
        if not isinstance(elt, O):
            self._count += 1

    @property
    def count(self):
        """Number of elements (including wires but not open-circuits)
        connected to the node"""

        return self._count

    def visible(self, draw_nodes):
        """Return true if node drawn"""

        if self.port:
            return True

        if draw_nodes in ('none', None, False):
            return False
        
        if draw_nodes == 'all':
            return True

        if draw_nodes == 'connections':
            return self.count > 2

        return self.name.find('_') == -1

    @property
    def port(self):
        """Return true if node is a port"""

        return self._port or self.count == 1


