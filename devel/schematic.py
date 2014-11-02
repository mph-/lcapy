from __future__ import print_function
import numpy as np


class Dnode(object):

    def __init__(self, name, port=False):

        self.name = name
        self.pos = None
        self.port = port

        parts = name.split('.')
        self.primary = len(parts) == 1


class NetElement(object):

    def __init__(self, name, node1, node2, symbol=None, orientation='up'):

        kind = name[0]
        if len(name) > 2 and name[0:2] == 'TF':
            kind = name[0:2]

        self.name = name
        self.symbol = symbol
        self.orientation = orientation
        self.nodes = (node1, node2)
        self.size = 2


    def __repr__(self):

        str = ', '.join(arg.__str__() for arg in [self.name] + list(self.nodes))
        return 'NetElement(%s)' % str


    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes])



class Schematic(object):

    def __init__(self, filename=None):

        self.elements = {}
        self.nodes = {}
        self.num_nodes = 0

        if filename is not None:
            self.netfile_add(filename)


    def __getitem__(self, name):
        """Return component by name"""

        return self.elements[name]


    def _nodeindex(self, node):
        """Return node index; ground is -1"""
        return self.revnodemap[node] - 1


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



    def _dnode_add(self, node, elt):

        if not self.nodes.has_key(node):
            self.nodes[node] = []
        self.nodes[node].append(elt)


    def _elt_add(self, elt):

        if self.elements.has_key(elt.name):
            print('Overriding component %s' % elt.name)     
            # Need to search lists and update component.
           
        self.elements[elt.name] = elt

        for dnode in elt.nodes:
            self._dnode_add(dnode, elt)
        

    def net_add(self, line):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        fields = line.split(';')

        parts = fields[0].split(' ')

        if len(fields) == 1:
            orientation = 'up'
        else:
            orientation = fields[1].strip()

        elt = NetElement(*parts, orientation=orientation)

        self._elt_add(elt)


    def _positions_calculate(self):

        num_nodes = len(self.nodes)

        A = np.zeros((num_nodes, num_nodes))
        bx = np.zeros(num_nodes)
        by = np.zeros(num_nodes)

        dnode_name_list = list(self.nodes)

        # Generate x and y constraint matrices and x and y component size vectors.
        k = 0
        for m, elt in enumerate(self.elements.values()):

            n1, n2 = elt.nodes[0], elt.nodes[1]
            m1, m2 = dnode_name_list.index(n1), dnode_name_list.index(n2)

            if k == 0:
                # Set first dnode to be arbitrary origin; this gets changed later.
                A[k, m1] = 1
                A[k, m1] = 1
                k += 1

            A[k, m1] = -1
            A[k, m2] = 1

            if elt.orientation == 'right':
                bx[k] = -elt.size
            elif elt.orientation == 'left':
                bx[k] = elt.size
            elif elt.orientation == 'up':
                by[k] = -elt.size
            elif elt.orientation == 'down':
                by[k] = elt.size
            else:
                raise ValueError('Unknown orientation %s' % elt.orientation)
            
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
#            print('%s @ (%.1f, %.1f)' % (dnode_name_list[m], x[m], y[m]))


        for m, elt in enumerate(self.elements.values()):

            n1, n2 = elt.nodes[0], elt.nodes[1]
            m1, m2 = dnode_name_list.index(n1), dnode_name_list.index(n2)

            elt.pos1 = pos[m1]
            elt.pos2 = pos[m2]

        self.dnode_positions = pos
        self.dnode_name_list = dnode_name_list


    def draw(self, filename=None):

        if not hasattr(self, 'dnode_positions'):
            self._positions_calculate()

        if filename != None:
            outfile = open(filename, 'w')
        else:
            import sys
            outfile = sys.stdout

        # Preamble
        print(r'\begin{tikzpicture}', file=outfile)

        # Write coordinates
        for m, dnode in enumerate(self.dnode_name_list):
            print(r'    \coordinate (%s) at (%.1f, %.1f);' % (dnode, self.dnode_positions[m][0], self.dnode_positions[m][1]), file=outfile)


        # Draw components
        for m, elt in enumerate(self.elements.values()):

            cpt = elt.name[0:1]

            # Need to special case port component.
            if cpt[0] == 'P':
                continue

            print(r'    \draw (%s) to [%s=$%s$] (%s);' % (elt.nodes[1], cpt, elt.name, elt.nodes[0]))

    
        # Label primary nodes
        if False:
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

    return sch

    
sch = test()
sch._positions_calculate()
sch.draw()
