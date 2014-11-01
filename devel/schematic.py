import numpy as np

def _node(name):

    parts = name.split('.')
    name = parts[0]
            
    # Convert to integer if possible
    try:
        node = int(name)
    except:
        node = name
    return node


class NetElement(object):

    def __init__(self, name, node1, node2, symbol=None, orientation='up'):

        kind = name[0]
        if len(name) > 2 and name[0:2] == 'TF':
            kind = name[0:2]

        self.name = name
        self.nodes = (_node(node1), _node(node2))
        self.symbol = symbol
        self.orientation = orientation
        self.dnodes = (node1, node2)


    def __repr__(self):

        str = ', '.join(arg.__str__() for arg in [self.name] + list(self.nodes))
        return 'NetElement(%s)' % str


    def __str__(self):

        return ' '.join(['%s' % arg for arg in (self.name, ) + self.nodes])



class Schematic(object):

    def __init__(self, filename=None):

        self.elements = {}
        self.nodes = {}
        self.dnodes = {}
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


    def _node_add(self, node, elt):

        if not self.nodes.has_key(node):
            self.nodes[node] = []
        self.nodes[node].append(elt)


    def _dnode_add(self, node, elt):

        if not self.dnodes.has_key(node):
            self.dnodes[node] = []
        self.dnodes[node].append(elt)


    def _elt_add(self, elt):

        if self.elements.has_key(elt.name):
            print('Overriding component %s' % elt.name)     
            # Need to search lists and update component.
           
        self.elements[elt.name] = elt

        self._node_add(elt.nodes[0], elt)
        self._node_add(elt.nodes[1], elt)

        self._dnode_add(elt.dnodes[0], elt)
        self._dnode_add(elt.dnodes[1], elt)
        

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


    def foo(self):

        num_dnodes = len(self.dnodes)

        A = np.zeros((num_dnodes, num_dnodes))
        #bx = np.zeros((num_dnodes, 1))
        #by = np.zeros((num_dnodes, 1))
        bx = np.zeros(num_dnodes)
        by = np.zeros(num_dnodes)

        dnodelist = list(self.dnodes)

        k = 0
        for m, elt in enumerate(self.elements.values()):
            print elt.name, elt.dnodes

            n1 = elt.dnodes[0]
            n2 = elt.dnodes[1]

            m1, m2 = dnodelist.index(n1), dnodelist.index(n2)

            if k == 0:
                # Set first dnode to be arbitrary origin; this gets changed later.
                A[k, m1] = 1
                A[k, m1] = 1
                k += 1

            A[k, m1] = -1
            A[k, m2] = 1

            if elt.orientation == 'right':
                bx[k] = -1
            elif elt.orientation == 'left':
                bx[k] = 1
            elif elt.orientation == 'up':
                by[k] = -1
            elif elt.orientation == 'down':
                by[k] = 1
            else:
                raise ValueError('Unknown orientation %s' % elt.orientation)
            
            print m1, m2
            k += 1

        Apinv = np.linalg.pinv(A)
        x = np.dot(Apinv, bx)
        y = np.dot(Apinv, by)

        x = x - x.min()
        y = y - y.min()
        
        print A
        #print bx
        #print by

        pos = np.zeros((num_dnodes, 2))
        for m in range(num_dnodes):
            pos[m][0] = x[m]
            pos[m][1] = y[m]
            print('%s @ (%.1f, %.1f)' % (dnodelist[m], x[m], y[m]))

        #print pos



def test():
    
    sch = Schematic()

    sch.net_add('P1 1 0.1')
    sch.net_add('R1 3 1; right')
    sch.net_add('L1 2 3; right')
    sch.net_add('C1 3 0; up')
    sch.net_add('P2 2 0.2')

    return sch

    
sch = test()
sch.foo()
