from __future__ import print_function
import numpy as np


class Node(object):

    def __init__(self, name):

        self.name = name
        self.pos = None
        self.port = False
        parts = name.split('_')
        self.rootname = parts[0]
        self.primary = len(parts) == 1
        self.list = []

    
    def append(self, elt):

        cpt = elt.name[0:1]        
        if cpt == 'P':
            self.port = True

        self.list.append(elt)


    @property
    def symbol(self):
        
        return 'o' if self.port else '*'



class NetElement(object):

    def __init__(self, name, node1, node2, symbol=None, dir='up', size=1):

        kind = name[0]
        if len(name) > 2 and name[0:2] == 'TF':
            kind = name[0:2]

        node1 = node1.replace('.', '_')
        node2 = node2.replace('.', '_')

        self.name = name
        self.symbol = symbol
        self.dir = dir
        self.nodes = (node1, node2)
        self.size = float(size) * 2


    def __repr__(self):

        str = ', '.join(arg.__str__() for arg in [self.name] + list(self.nodes))
        return 'NetElement(%s)' % str


    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes])



class Schematic(object):

    def __init__(self, filename=None):

        self.elements = {}
        self.nodes = {}
        self.vnodes = {}

        if filename is not None:
            self.netfile_add(filename)


    def __getitem__(self, name):
        """Return component by name"""

        return self.elements[name]


    def netfile_add(self, filename):    
        """Add the nets from file with specified filename"""

        file = open(filename, 'r')
        
        lines = file.readlines()

        for line in lines:
            # Skip comments
            if line[0] in ('#', '%'):
                continue
            self.net_add(line.strip())


    def netlist(self):
        """Return the current netlist"""

        return '\n'.join([elt.__str__() for elt in self.elements.values()])


    def _node_add(self, node, elt):

        if not self.nodes.has_key(node):
            self.nodes[node] = Node(node)
        self.nodes[node].append(elt)

        vnode = self.nodes[node].rootname
        if not self.vnodes.has_key(vnode):
            self.vnodes[vnode] = {}
        if not self.vnodes[vnode].has_key(node):
            self.vnodes[vnode][node] = elt


    def _elt_add(self, elt):

        if self.elements.has_key(elt.name):
            print('Overriding component %s' % elt.name)     
            # Need to search lists and update component.
           
        self.elements[elt.name] = elt

        for node in elt.nodes:
            self._node_add(node, elt)
        

    def _hints_parse(self, str):

        hints = {'dir' : 'up',
                 'size' : 1}

        for part in str.split(','):
            part = part.strip()

            if part in ('up', 'down', 'left', 'right'):
                hints['dir'] = part
                continue

            fields = part.split('=')
            key = fields[0].strip()
            arg = fields[1].strip() if len(fields) > 1 else ''
            hints[key] = arg

        return hints


    def net_add(self, line):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        fields = line.split(';')

        str = fields[1] if len(fields) > 1 else ''

        hints = self._hints_parse(str)

        parts = fields[0].split(' ')
        elt = NetElement(*parts, dir=hints['dir'], size=hints['size'])

        self._elt_add(elt)


    def _positions_calculate(self):

        num_nodes = len(self.nodes)

        A = np.zeros((num_nodes, num_nodes))
        bx = np.zeros(num_nodes)
        by = np.zeros(num_nodes)

        node_name_list = list(self.nodes)

        # Generate x and y constraint matrices and x and y component size vectors.
        k = 0
        for m, elt in enumerate(self.elements.values()):

            n1, n2 = elt.nodes[0], elt.nodes[1]
            m1, m2 = node_name_list.index(n1), node_name_list.index(n2)

            if k == 0:
                # Set first node to be arbitrary origin; this gets changed later.
                A[k, m1] = 1
                A[k, m1] = 1
                k += 1

            A[k, m1] = -1
            A[k, m2] = 1

            if elt.dir == 'right':
                bx[k] = -elt.size
            elif elt.dir == 'left':
                bx[k] = elt.size
            elif elt.dir == 'up':
                by[k] = -elt.size
            elif elt.dir == 'down':
                by[k] = elt.size
            else:
                raise ValueError('Unknown dir %s' % elt.dir)
            
            k += 1

        Apinv = np.linalg.pinv(A)
        x = np.dot(Apinv, bx)
        y = np.dot(Apinv, by)

        # Adjust positions so origin at (0, 0).
        x = x - x.min()
        y = y - y.min()
        
        #print A
        #print bx
        #print by

        pos = np.zeros((num_nodes, 2))
        for m in range(num_nodes):
            pos[m][0] = x[m]
            pos[m][1] = y[m]
#            print('%s @ (%.1f, %.1f)' % (node_name_list[m], x[m], y[m]))


        for m, elt in enumerate(self.elements.values()):

            n1, n2 = elt.nodes[0], elt.nodes[1]
            m1, m2 = node_name_list.index(n1), node_name_list.index(n2)

            elt.pos1 = pos[m1]
            elt.pos2 = pos[m2]

        self.node_positions = pos
        self.node_name_list = node_name_list


    def draw(self, draw_nodes=True, label_nodes=True, filename=None):

        if not hasattr(self, 'node_positions'):
            self._positions_calculate()

        if filename != None:
            outfile = open(filename, 'w')
        else:
            import sys
            outfile = sys.stdout

        # Preamble
        print(r'\begin{tikzpicture}', file=outfile)

        # Write coordinates
        for m, node in enumerate(self.node_name_list):
            print(r'    \coordinate (%s) at (%.1f, %.1f);' % (node, self.node_positions[m][0], self.node_positions[m][1]), file=outfile)


        # Draw components
        for m, elt in enumerate(self.elements.values()):

            cpt = elt.name[0:1]

            # Need to special case port component.
            if cpt[0] == 'P':
                continue

            node_str = ''
            if draw_nodes:
                node_str = self.nodes[elt.nodes[1]].symbol + '-' + self.nodes[elt.nodes[0]].symbol

            print(r'    \draw (%s) to [%s=$%s$, %s] (%s);' % (elt.nodes[1], cpt, elt.name, node_str, elt.nodes[0]))


        # Draw wires
        for m, vnode in enumerate(self.vnodes.values()):
            num_wires = len(vnode) - 1
            if num_wires == 0:
                continue
            # TODO: remove overdraw wires...

            for n in range(num_wires):
                n1 = vnode.keys()[n]
                n2 = vnode.keys()[n + 1]

                node_str = ''
                if draw_nodes:
                    node_str = self.nodes[n1].symbol + '-' + self.nodes[n2].symbol

                    print(r'    \draw (%s) to [short, %s] (%s);' % (n1, node_str, n2))
    
        # Label primary nodes
        if label_nodes:
            for m, node in enumerate(self.nodes.values()):
                if not node.primary:
                    continue
                print(r'    \draw {[anchor=south east] (%s) node {%s}};' % (node.name, node.name))

        print(r'\end{tikzpicture}', file=outfile)


def test():
    
    sch = Schematic()

    sch.net_add('P1 1 0.1')
    sch.net_add('R1 3 1; right')
    sch.net_add('L1 2 3; right')
    sch.net_add('C1 3 0; up')
    sch.net_add('P2 2 0.2')

    sch.draw()
    return sch
