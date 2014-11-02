from __future__ import print_function
import numpy as np
import re


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

    def __init__(self, name, node1, node2, symbol=None, opts=None):

        # Regular expression alternate matches stop with first match
        # so need to have longer ones first.
        cpts = ('R', 'C', 'L', 'Vac', 'Vdc', 'Iac', 'Idc', 'V', 'I',
                'TF', 'P', 'port', 'W', 'wire')


        match = re.match(r'(%s)(\w)' % '|'.join(cpts), name)
        if not match:
            raise ValueError('Unknown component %s' % name)

        cpt = match.groups(1)[0]
        id = match.groups(1)[1]

        node1 = node1.replace('.', '_')
        node2 = node2.replace('.', '_')

        self.symbol = symbol

        if symbol is None:
            symbol = cpt + '_{' + id + '}'

        self.name = name
        self.cpt = cpt
        self.autosymbol = symbol
        self.nodes = (node1, node2)
        self.opts = opts


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
        self.scale = 2

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
        

    def _opts_parse(self, str):

        opts = {'dir' : 'up',
                 'size' : 1}

        for part in str.split(','):
            part = part.strip()

            if part in ('up', 'down', 'left', 'right'):
                opts['dir'] = part
                continue

            fields = part.split('=')
            key = fields[0].strip()
            arg = fields[1].strip() if len(fields) > 1 else ''
            opts[key] = arg

        return opts


    def net_add(self, line):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive nose and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        fields = line.split(';')

        str = fields[1] if len(fields) > 1 else ''

        opts = self._opts_parse(str)

        parts = fields[0].split(' ')
        elt = NetElement(*parts, opts=opts)

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

            dir = elt.opts['dir']
            size = float(elt.opts['size']) * self.scale

            if dir == 'right':
                bx[k] = -size
            elif dir == 'left':
                bx[k] = size
            elif dir == 'up':
                by[k] = -size
            elif dir == 'down':
                by[k] = size
            else:
                raise ValueError('Unknown dir %s' % dir)
            
            k += 1

        Apinv = np.linalg.pinv(A)
        x = np.dot(Apinv, bx)
        y = np.dot(Apinv, by)

        # Adjust positions so origin at (0, 0).
        x = x - x.min()
        y = y - y.min()

        self.xcentre = x.mean()
        self.ycentre = y.mean()
        
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


    def _make_wires1(self, vnode):

        num_wires = len(vnode) - 1
        if num_wires == 0:
            return []

        wires = []

        # TODO: remove overdraw wires...
        for n in range(num_wires):
            n1 = vnode.keys()[n]
            n2 = vnode.keys()[n + 1]
            
            wires.append(NetElement('W_', n1, n2))

        return wires


    def _make_wires(self):

        wires = []

        vnode_dir = self.vnodes

        for m, vnode in enumerate(vnode_dir.values()):
            wires.extend(self._make_wires1(vnode))
            
        return wires


    def tikz_draw(self, draw_labels=True, draw_nodes=True, label_nodes=True, filename=None, args=None):

        self._positions_calculate()

        if filename != None:
            outfile = open(filename, 'w')
        else:
            import sys
            outfile = sys.stdout

        # Preamble
        if args is None: args = ''
        print(r'\begin{tikzpicture}[%s]' % args, file=outfile)

        # Write coordinates
        for m, node in enumerate(self.node_name_list):
            print(r'    \coordinate (%s) at (%.1f, %.1f);' % (node, self.node_positions[m][0], self.node_positions[m][1]), file=outfile)


        cpt_map = {'R' : 'R', 'C' : 'C', 'L' : 'L', 'V' : 'V', 'I' : 'I',
                   'Vac' : 'sV', 'Vdc' : 'V', 'Iac' : 'sI', 'Idc' : 'I', 
                   'TF' : 'transformer', 'P' : 'open', 'port' : 'open',
                   'W' : 'short', 'wire' : 'short'}

        # Draw components
        for m, elt in enumerate(self.elements.values()):

            n1 = elt.nodes[1]
            n2 = elt.nodes[0]
            cpt = cpt_map[elt.cpt]

            # Need to special case port component.
            if cpt[0] == 'P':
                if elt.symbol:
                    dir = '^' if elt.pos1[0] > self.xcentre else ''

                    print(r'    \draw (%s) to [open, v%s=$%s$] (%s);' % (n2, dir, elt.symbol, n1))
                else:
                    print(r'    \draw (%s) to [open] (%s);' % (n2, n1))
                continue


            # Current, voltage, label options.
            # It might be better to allow any options and prune out
            # dir and size.
            opts_str = ''
            for opt in ('i', 'i_', 'i^', 'i_>', 'i_<', 'i^>', 'i^<', 'v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<', 'l', 'l^', 'l_'):
                if elt.opts.has_key(opt):
                    opts_str += '%s=$%s$, ' % (opt, elt.opts[opt])

            node_str = ''
            if draw_nodes:
                node_str = self.nodes[n1].symbol + '-' + self.nodes[n2].symbol

               
            label_str =''
            if draw_labels and not ('l' in elt.opts.keys() or 'l_' in elt.opts.keys() or 'l^' in elt.opts.keys()):
                label_str = '=$%s$' % elt.autosymbol
            

            print(r'    \draw (%s) to [%s%s, %s%s] (%s);' % (n1, cpt, label_str, opts_str, node_str, n2))

        wires = self._make_wires()

        # Draw wires
        for wire in wires:
            n1 = wire.nodes[1]
            n2 = wire.nodes[0]

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


    def draw(self, draw_labels=True, draw_nodes=True, label_nodes=True,
             filename=None, args=None, scale=2):

        self.scale = scale

        return self.tikz_draw(draw_labels=draw_labels, draw_nodes=draw_nodes,
                              label_nodes=label_nodes, filename=filename, args=args)



def test():
    
    sch = Schematic()

    sch.net_add('P1 1 0.1')
    sch.net_add('R1 3 1; right')
    sch.net_add('L1 2 3; right')
    sch.net_add('C1 3 0; up')
    sch.net_add('P2 2 0.2')

    sch.draw()
    return sch
