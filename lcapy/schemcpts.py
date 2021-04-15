"""
This module defines and draws the schematic components using
circuitikz.   The components are defined at the bottom of this file.

Copyright 2015--2021 Michael Hayes, UCECE
"""


from __future__ import print_function
from .latex import latex_format_label
from .schemmisc import Pos, Steps
from .opts import Opts
import numpy as np
import sys

module = sys.modules[__name__]

# There are two types of component (Cpt).
#
# 1. The stretchable components (such as resistors and capacitors) have
# wires that can be stretched.  The size attribute controls the
# spacing between the nodes but does not affect the component size.
# The component size can be changed with the scale attribute (this changes
# the dipole length).  The aspect ratio is not easy to change (need to use
# dipole/resistor/height).
#
# 2.  The fixed components (such as ICs) do not have wires and cannot
# be stretched.  The (unrotated) width is set by the size attribute and
# the (unrotated) height is set by the aspect attribute.  The scale attribute
# is not used.

# There are two paradigms used for specifying node coordinates:
#
# 1.  The old model.  The required_node_names method returns a subset of
# node_names as a list.
#
# 2.  The new model.  The node_pinnames attribute specifies the subset
# of required nodes.  It is a list of pinnames that are used to find
# the pin coordinates and thus the node coordinates.  Note, some
# components to do not have any explicit nodes (shapes, chips, etc).

# The direction commands rotate the component:
# right (0 degrees)
# up    (90 degrees)
# left  (180 degrees)
# down  (-90 degrees)
#
# Lcapy uses this information to position the second node with respect
# to the first.   Circuitikz then infers the rotation.

def check_boolean(value):

    if value not in (True, False, 'none', None, 'true', 'false'):
        raise ValueError('Unexpected Boolean value %s' % value)
    return value in (True, 'true')


class Cpt(object):

    voltage_keys = ('v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<',
                    'v<', 'v>')
    current_keys = ('i', 'i_', 'i^', 'i_>',  'i_<', 'i^>', 'i^<',
                    'i>_', 'i<_', 'i>^', 'i<^', 'i>', 'i<', 'ir')
    flow_keys = ('f', 'f_', 'f^', 'f_>',  'f_<', 'f^>', 'f^<',
                    'f>_', 'f<_', 'f>^', 'f<^', 'f>', 'f<')    
    label_keys = ('l', 'l_', 'l^')
    inner_label_keys = ('t', )    
    implicit_keys =  ('implicit', 'ground', 'sground', 'rground',
                      'cground', 'nground', 'pground', 'vss', 'vdd',
                      'vee', 'vcc', 'input', 'output', 'bidir', 'pad',
                      'antenna', 'rxantenna', 'txantenna')
    # The following keys do not get passed through to circuitikz.
    misc_keys = ('left', 'right', 'up', 'down', 'rotate', 'size',
                 'mirror', 'invert', 'scale', 'invisible', 'variable', 'fixed',
                 'aspect', 'pins', 'image', 'offset', 'pinlabels',
                 'pinnames', 'pinnodes', 'pindefs', 'outside',
                 'pinmap', 'kind', 'wire', 'ignore', 'style',
                 'nowires', 'nolabels', 'steps', 'free', 'fliplr', 'flipud',
                 'nodots', 'draw_nodes', 'label_nodes', 'nodraw',
                 'mirrorinputs')

    can_rotate = True
    can_scale = False
    can_mirror = False
    can_invert = False
    do_transpose = False    
    can_stretch = True
    default_width = 1.0
    default_aspect = 1.0
    # node_pinnames maps node numbers to pinnames
    node_pinnames = ()
    default_pins = None
    pins = {}
    # Auxiliary nodes are used for finding the centre of the shape or
    # to define a bounding box.
    auxiliary = {}
    directive = False
    place = True

    @property
    def s(self):
        """Sanitised name"""
        return self.name.replace('.', '@')

    def __init__(self, sch, namespace, defname, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        self.sch = sch
        self.type = cpt_type
        self.id = cpt_id
        self.defname = defname        
        self.name = name
        self.namespace = namespace

        self.string = string
        self.net = string.split(';')[0]
        self.opts_string = opts_string
        self.opts = Opts(self.opts_string)        

        self.args = args
        self.classname = self.__class__.__name__

        # Drawing hints
        self.opts = Opts(opts_string)

        # The ordering of this list is important.
        self.node_names = node_names

        # The relative node names start with a .
        self.relative_node_names = []
        for name in node_names:
            fields = name.split('.')
            if len(fields) < 2:
                continue
            if fields[-2] == self.name:
                self.relative_node_names.append(name)

        prefix = self.name + '.'
        
        auxiliary_node_names = []
        for pin in self.auxiliary:
            auxiliary_node_names.append(prefix + pin)

        self.auxiliary_node_names = auxiliary_node_names

        self.allpins = self.pins.copy()
        self.allpins.update(self.auxiliary)

        pin_node_names = []
        for pin in self.pins.keys():
            pin_node_names.append(prefix + pin)
        
        self.pindefs = self.parse_pindefs()
        for nodename, pindef in self.pindefs.items():
            if True:
                # Remove old name.
                self.allpins[pindef] = self.allpins.pop(nodename)
                pin_node_names.remove(prefix + nodename)
            else:
                # Keep old name as an alias.  This will cause problems
                # when printing pinnodes.
                self.allpins[pindef] = self.allpins[nodename]                
            pin_node_names.append(prefix + pindef)            
        
        # These are all the pin names belonging to the cpt.
        self.pin_node_names = pin_node_names
            
        # These are all the pin nodes required to be shown for the cpt.
        # This is set by the process_pins method.
        self.drawn_pins = []

        self.all_node_names = list(self.required_node_names) + auxiliary_node_names + pin_node_names

        # Create dictionary of pinnames sharing the same relative
        # coords (pinname aliases).
        coords = {}
        for pinname, coord in self.pins.items():
            if coord not in coords:
                coords[coord] = [pinname]
            else:
                coords[coord] += [pinname]

        self.pinname_coords = coords
        if False:
            for coord, pinnames in coords.items():
                if len(pinnames) > 1:
                    print('%s: pinname aliases: %s' % (self.name, pinnames))
        
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.string

    @property
    def size(self):
        """Component size between its nodes"""
        if 'size' in self.opts:
            val = self.opts['size']
        elif self.right:
            val = self.opts['right']
        elif self.down:
            val = self.opts['down']
        elif self.left:
            val = self.opts['left']
        elif self.up:
            val = self.opts['up']
        else:
            val = self.default_width
        if val == '':
            val = self.default_width
        return float(val)

    @property
    def scale(self):
        return float(self.opts.get('scale', 1.0))

    @property
    def down(self):
        return 'down' in self.opts

    @property
    def up(self):
        return 'up' in self.opts

    @property
    def left(self):
        return 'left' in self.opts

    @property
    def right(self):
        return 'right' in self.opts

    @property
    def horizontal(self):
        return self.left or self.right

    @property
    def vertical(self):
        return self.up or self.down

    def boolattr(self, opt):

        if opt not in self.opts:
            return False
        if self.opts[opt] == '':
            return True
        return self.opts[opt]

    @property
    def mirrorinputs(self):
        return self.boolattr('mirrorinputs')
    
    @property
    def mirror(self):
        return self.boolattr('mirror') or self.flipud

    @property
    def invert(self):
        return self.boolattr('invert')  or self.fliplr

    @property
    def flipud(self):
        return self.boolattr('flipud')

    @property
    def fliplr(self):
        return self.boolattr('fliplr')    

    @property
    def nodots(self):
        return self.boolattr('nodots')

    @property
    def nolabels(self):
        return self.boolattr('nolabels')
   
    @property
    def wire(self):
        return self.boolattr('wire')

    @property
    def nowires(self):
        return self.boolattr('nowires')        

    @property
    def invisible(self):
        return self.boolattr('invisible')

    @property
    def nodraw(self):
        return self.boolattr('nodraw')    

    @property
    def ignore(self):
        return self.boolattr('ignore')    

    @property
    def variable(self):
        return self.boolattr('variable')

    @property
    def fixed(self):
        return self.boolattr('fixed')

    @property
    def stretch(self):
        return self.can_stretch and not self.fixed

    @property
    def angle(self):
        """Return rotation angle"""
        if self.right:
            angle = 0   
        elif self.down:
            angle = -90
        elif self.left:
            angle = 180
        elif self.up:
            angle = 90
        else:
            angle = -90 if self.type in ('P', ) else 0
        
        if 'rotate' in self.opts:
            angle += float(self.opts['rotate'])
        return angle

    @property
    def w(self):
        """Normalised width"""
        return 1.0

    @property
    def h(self):
        """Normalised height"""
        return self.w / self.aspect

    @property
    def aspect(self):
        return float(self.opts.get('aspect', self.default_aspect))

    @property
    def offset(self):
        return float(self.opts.get('offset', 0))

    @property
    def steps(self):
        return self.opts.get('steps', None)

    @property
    def free(self):
        return self.boolattr('free')

    @property
    def kind(self):
        return self.opts.get('kind', None)

    @property
    def draw_nodes_opt(self):
        return self.opts.get('draw_nodes', None)

    @property
    def label_nodes_opt(self):
        return self.opts.get('label_nodes', None)    

    @property
    def style(self):
        return self.opts.get('style', None)
    
    def R(self, angle_offset=0):
        """Return rotation matrix"""
        angle = self.angle + angle_offset
        
        Rdict = {0: ((1, 0), (0, 1)),
                 90: ((0, 1), (-1, 0)),
                 180: ((-1, 0), (0, -1)),
                 -180: ((-1, 0), (0, -1)),
                 -90: ((0, -1), (1, 0))}
        if angle in Rdict:
            return np.array(Rdict[angle])
        
        t = angle / 180.0 * np.pi
        return np.array(((np.cos(t), np.sin(t)),
                         (-np.sin(t), np.cos(t))))

    @property
    def required_node_names(self):
        """Subset of node_names.  This filters out nodes that are not
        drawn.  For example, the ground node of an Eopamp is not drawn."""

        # Old model.   The number of node names in the list
        # must match the number of entries in coords.
        if self.node_pinnames == ():
            return self.node_names

        # New model.  The node_pinnames tuple specifies the required
        # nodes.
        node_names = []
        for pinname, node_name in zip(self.node_pinnames, self.node_names):
            if pinname != '':
                node_names.append(node_name)
        return tuple(node_names)

    def node(self, pinname):
        """Return node by pinname"""

        if pinname in self.node_pinnames:
            index = self.node_pinnames.index(pinname)
            node_name = self.node_names[index]
        else:
            node_name = self.name + '.' + pinname
        for node in self.nodes:
            if node.name == node_name:
                return node
        raise ValueError('Unknown pinname %s for %s' % (pinname, self))

    def pinpos(self, pinname):
        """Return pinpos by pinname"""

        if pinname not in self.allpins:
            raise ValueError('Unknown pin %s for %s, known pins: %s'
                             % (pinname, self.name, ', '.join(list(self.allpins))))
        return self.allpins[pinname][0][0]

    def fakepin(self, pinname):
        """Return True if pinname is a fake pin"""

        # Check if known pinname; it might be a netlist node name
        if pinname not in self.allpins:
            return True

        # Fake pins have a pinpos like lx, rx, tx, or bx
        return len(self.allpins[pinname][0]) == 2
    
    @property
    def nodes(self):
        """Nodes used to draw the current element."""

        if hasattr(self, '_nodes'):
            return self._nodes

        # Perhaps determine coords here as well and cache them?

        # The ordering of the first nodes is important.   
        # FIXME: there can be duplicates but cannot use a set to
        # remove them since this will change the ordering.
        node_names = self.all_node_names
        
        rnodes = []
        for n in node_names:
            if n in self.sch.nodes:
                rnodes.append(self.sch.nodes[n])
        self._nodes = rnodes
        return rnodes

    @property
    def drawn_nodes(self):
        """Nodes that are drawn"""
        return self.nodes

    @property
    def required_pins(self):

        rpins = []
        for node in self.nodes:
            node_name = node.name
            if node_name in self.node_names:
                index = self.node_names.index(node_name)
                pinname = self.node_pinnames[index]
            elif node_name in self.sch.nodes:
                pinname = node_name.split('.')[-1]
            else:
                raise ValueError('Unknown node %s' % node_name)

            if pinname != '':
                rpins.append(self.allpins[pinname])
        return rpins

    @property
    def coords(self):

        rpins = self.required_pins
        coords = [pin[1:] for pin in rpins]
        return coords

    @property
    def scales(self):

        if not self.can_scale:
            return [1] * len(self.coords)

        rpins = self.required_pins
        scales = []
        for pin in rpins:
            scale = 1
            pinpos = pin[0]
            if not pinpos.endswith('x'):
                scale = 2 * self.scale / self.width
            scales.append(scale)
        return scales

    @property
    def tcoords(self):
        """Transformed coordinates for each of the nodes"""

        if hasattr(self, '_tcoords'):
            return self._tcoords

        coords = self.coords
        scales = self.scales

        tcoords = np.zeros((len(coords), 2))
        for m in range(len(coords)):
            tcoords[m] = self.tf((0, 0), coords[m], scale=scales[m])
        self._tcoords = tcoords
        return self._tcoords

    @property
    def xvals(self):
        return self.tcoords[:, 0]

    @property
    def yvals(self):
        return self.tcoords[:, 1]

    def midpoint(self, node1, node2):
        return (node1.pos + node2.pos) * 0.5

    def _node_str(self, node1, **kwargs):

        draw_nodes = self.draw_nodes_opt
        if draw_nodes is None:
            draw_nodes = kwargs.get('draw_nodes', True)

        s = ''
        if node1.visible(draw_nodes) and not node1.pin:
            s = 'o' if node1.port else '*'
        return s

    def _node_pair_str(self, node1, node2, **kwargs):

        # Create o-o o-* *-* etc.
        s = self._node_str(node1, **kwargs)
        s += '-'
        s += self._node_str(node2, **kwargs)

        if s == '-':
            s = ''
        
        return s

    def draw_node(self, n, **kwargs):

        # Don't draw nodes for open-circuits.  Use port
        # if want nodes drawn.
        if self.type == 'O':
            return ''
        
        draw_nodes = self.draw_nodes_opt
        if draw_nodes is None:
            draw_nodes = kwargs.get('draw_nodes', True)

        s = ''
        if not draw_nodes:
            return s

        if not n.visible(draw_nodes) or n.pin:
            return s

        if n.port:
            s = r'  \draw (%s) node[ocirc] {};''\n' % n.s
        else:
            s = r'  \draw (%s) node[circ] {};''\n' % n.s

        return s

    def draw_nodes(self, **kwargs):

        s = ''
        for n in self.drawn_nodes:
            s += self.draw_node(n, **kwargs)
        return s

    def draw_pins(self):

        s = ''
        for n in self.drawn_pins:        
            s += r'  \draw (%s) node[ocirc] {};''\n' % n.s
        return s

    def draw_pinlabel(self, node):

        if node.pinlabel == '':
            return ''

        pinpos = node.pinpos
        pinpos = self.pinpos_rotate(pinpos, self.angle)
        
        outside = self.opts.get('outside', False)

        if self.invert:
            if pinpos == 'l':
                pinpos = 'r'
            elif pinpos == 'r':
                pinpos = 'l'
        if self.mirror:
            if pinpos == 't':
                pinpos = 'b'
            elif pinpos == 'b':
                pinpos = 't'                                
        
        if outside == '':            
            anchors = {None: 'south east',
                       'c': 'south east',
                       'l' : 'south east', 'r' : 'north west', 
                       't' : 'south west', 'b' : 'north west'}
        else:
            anchors = {None: 'south east',
                       'c': 'south east',
                       'l' : 'west', 'r' : 'east', 
                       't' : 'north', 'b' : 'south'}

        anchor = anchors[pinpos]

        return r'  \draw[anchor=%s] (%s) node {%s};''\n' % (
            anchor, node.s, node.pinlabel.replace('_', r'\_'))

    def draw_pinname(self, node):

        if node.pinname == '':
            return ''

        pinpos = node.pinpos        
        pinpos = self.pinpos_rotate(pinpos, self.angle)

        # Move pinnames to outside
        mapping = {None: None,
                   'c': 'c',
                   'l' : 'r', 'r' : 'l', 
                   't' : 'b', 'b' : 't'}
        pinpos = mapping[pinpos]
            
        anchors = {None: 'south east',
                   'c': 'south east',
                   'l' : 'west', 'r' : 'east', 
                   't' : 'north', 'b' : 'south'}
        anchor = anchors[pinpos]

        return r'  \draw[anchor=%s] (%s) node {%s};''\n' % (
            anchor, node.s, node.pinname.replace('_', r'\_'))        
    
    def draw_node_label(self, node, label_nodes):

        if not node.show_label(label_nodes):
            return ''
            
        anchors = {None: 'south east',
                   'c': 'south east',                   
                   'l' : 'west', 'r' : 'east', 
                   't' : 'north', 'b' : 'south'}
        anchor = anchors[node.labelpos]

        return r'  \draw[anchor=%s] (%s) node {%s};''\n' % (
            anchor, node.s, node.label)
       
    def draw_node_labels(self, **kwargs):

        label_nodes = self.label_nodes_opt
        if label_nodes is None:
            label_nodes = kwargs.get('label_nodes', 'primary')

        s = ''
        for node in self.drawn_nodes:

            if node.pin:
                if not node.belongs(self.name):
                    continue
                s += self.draw_pinname(node)
                s += self.draw_pinlabel(node)

            else:
                if node.auxiliary:
                    continue
                s += self.draw_node_label(node, label_nodes)                

        return s
    
    def draw(self, **kwargs):
        raise NotImplementedError('draw method not implemented for %s' % self)

    def find_ref_node_names(self):
        """Determine which nodes are referenced for this component."""

        ref_node_names = []

        for nodename, node in self.sch.nodes.items():
            if nodename in self.relative_node_names:
                node.ref = self
            elif node.belongs(self.name):
                if node.basename not in self.allpins:
                    continue
                node.ref = self
                node.pin = not self.fakepin(node.basename)
                if node.basename not in self.auxiliary:
                    ref_node_names.append(node.name)
            elif self.namespace != '' and nodename.startswith(self.namespace):
                # Need to be lenient here since can have any old name.
                node.ref = self                

        return ref_node_names
    
    def setup(self):
        self.ref_node_names = self.find_ref_node_names()

    def parse_pindefs(self):
        return {}
        
    def opts_str_list(self, choices):
        """Format voltage, current, or label string as a key-value pair
        and return list of strings"""

        def fmt(key, val):
            label = latex_format_label(val)
            if label == '':
                label = '{}'
            if not (label[0] == '{' and label[-1] == '}'):
                label = '{' + label + '}'

            return '%s=%s' % (key, label)

        return [fmt(key, val) for key, val in self.opts.items()
                if key in choices]

    def opts_str(self, choices):
        """Format voltage, current, or label string as a key-value pair"""

        return ','.join(self.opts_str_list(choices))

    @property
    def voltage_str(self):

        return self.opts_str(self.voltage_keys)

    @property
    def current_str(self):

        if 'ir' in self.opts:
            label = self.opts.pop('ir')
            self.opts.add('i_<=' + label)
        
        return self.opts_str(self.current_keys)

    @property
    def flow_str(self):

        return self.opts_str(self.flow_keys)    

    @property
    def label_str(self):

        return self.opts_str(self.label_keys)

    @property
    def label_str_list(self):

        return self.opts_str_list(self.label_keys)

    @property
    def inner_label_str(self):

        return self.opts_str(self.inner_label_keys)    

    @property
    def args_list(self):

        def fmt(key, val):
            return '%s=%s' % (key, latex_format_label(val))

        keys = self.voltage_keys + self.current_keys + self.flow_keys + self.label_keys + self.inner_label_keys + self.misc_keys + self.implicit_keys
    
        return [fmt(key, val) for key, val in self.opts.items() if key not in keys]

    @property
    def args_str(self):

        return ','.join(self.args_list)
    
    def label(self, keys=('l', 'l^', 'l_'), default=True, **kwargs):

        label_values = kwargs.get('label_values', True)
        label_values = check_boolean(label_values)
        label_ids = kwargs.get('label_ids', True)
        label_ids = check_boolean(label_ids)        

        label_str = ''
        if label_ids is True:
            label_str = self.id_label
        if label_values and self.value_label != '':
            label_str = self.value_label        

        has_label = False
        for key, val in self.opts.items():
            if key in ('l', 'l^', 'l_'):
                has_label = True
                break

        if has_label:
            # Override label if specified.  There are no placement options.
            label_str = ','.join([latex_format_label(val)
                                  for key, val in self.opts.items()
                                  if key in keys])
        elif not default:
            return ''

        # Remove curly braces.
        if len(label_str) > 1 and label_str[0] == '{' and label_str[-1] == '}':
            label_str = label_str[1:-1]
        return label_str

    def label_make(self, label_pos='', **kwargs):

        # TODO merge with label

        label_values = kwargs.get('label_values', True)
        label_values = check_boolean(label_values)                
        label_ids = kwargs.get('label_ids', True)
        label_ids = check_boolean(label_ids)        

        # Generate default label.
        if (label_ids and label_values and self.id_label != '' 
            and self.value_label and self.id_label != self.value_label):
            label_str = r'l%s={%s}{=%s}' % (label_pos, self.id_label,
                                            self.value_label)
        elif label_ids and self.id_label != '':
            label_str = r'l%s=%s' % (label_pos, self.id_label)
        elif label_values and self.value_label != '':
            label_str = r'l%s=%s' % (label_pos, self.value_label)
        else:
            label_str = ''

        # Override label if specified.
        if self.label_str != '':
            label_str = self.label_str

        if self.inner_label_str != '':
            label_str += ',t=' + self.inner_label_str
        return label_str

    def check(self):
        """Check schematic options and return True if component is to be drawn"""

        if not self.can_rotate and self.angle != 0:
            raise ValueError('Cannot rotate component %s' % self.name)

        if not self.can_scale and self.scale != 1:
            raise ValueError('Cannot scale component %s' % self.name)

        if not self.can_mirror and (self.mirror or self.mirrorinputs):
            raise ValueError('Cannot mirror component %s' % self.name)

        if not self.can_invert and self.invert:
            raise ValueError('Cannot invert component %s' % self.name)        

        if self.left + self.right + self.up + self.down > 1:
            raise ValueError('Mutually exclusive drawing directions for %s' % self.name)

        return not self.invisible

    def tf(self, centre, offset, angle_offset=0.0, scale=None):
        """Transform coordinate"""

        raise NotImplementedError('tf method not implemented for %s' % self)

    def draw_path(self, points, style='', join='--', closed=False):

        path = (' %s ' % join).join(['(%s)' % point for point in points])
        if closed:
            path += ' %s cycle' % join

        args_str = self.args_str
        if style == '':
            s = args_str
        elif args_str == '':
            s = style
        else:
            s = style + ', ' + args_str
        if s != '':
            s = '[%s]' % s

        return r'  \draw%s %s;''\n' % (s, path)

    def annotate(self, pos, label, args_str='', bold=False):

        if bold:
            if label.startswith('$') and label.endswith('$'):
                # boldsymbol but requires amssym or amsmath.
                label = r'$\Large \boldsymbol{%s}$' % label[1:-1]
            else:
                label = r'\textbf{%s}' % label
        
        return r'  \draw[%s] (%s) node[] {%s};''\n'% (
            args_str, pos, label)        
    
    def draw_label(self, pos, keys=('l', 'l^', 'l_'), default=True, **kwargs):

        return self.annotate(pos, self.label(keys, default=default,
                                             **kwargs), self.args_str)
    

class A(Cpt):
    """Annotation."""

    place = False

    def draw(self, **kwargs):

        n = self.nodes[0]
        s = r'  \draw[%s] (%s) node {%s};''\n' % (self.args_str, n.s,
                                                  self.label(**kwargs))
        return s
    
        
class StretchyCpt(Cpt):

    can_stretch = True

    def xtf(self, centre, offset, angle_offset=0.0):
        """Transform coordinate."""

        # Note the size attribute is not used.   This only scales the x coords.
        if isinstance(offset[0], tuple):
            return [self.xtf(centre, offset1, angle_offset) for offset1 in offset]

        return centre + np.dot((offset[0] * self.w * self.scale, offset[1] * self.h), self.R(angle_offset)) * self.sch.node_spacing


    def tf(self, centre, offset, angle_offset=0.0, scale=None):
        """Transform coordinate."""

        # Note the size attribute is not used.
        if hasattr(offset[0], '__iter__'):        
            return [self.tf(centre, offset1, angle_offset, scale) for offset1 in offset]
        x, y = offset

        if self.do_transpose:
            if self.mirror:
                y = -y
            if self.invert:
                x = -x
        
        if scale is None:
            scale = self.scale * self.sch.node_spacing
        
        return centre + np.dot((x * self.w, y * self.h), self.R(angle_offset)) * scale


class FixedCpt(Cpt):

    can_stretch = False

    @property
    def centre(self):
        if hasattr(self, 'pins'):
            # Look for centre pin.
            for node in self.nodes:
                if node.name.split('.')[-1] == 'mid':
                    return node.pos
        
        N = len(self.nodes)
        return self.midpoint(self.nodes[0], self.nodes[N // 2])

    def tf(self, centre, offset, angle_offset=0.0, scale=None):
        """Transform coordinate."""

        if hasattr(offset[0], '__iter__'):
            return [self.tf(centre, offset1, angle_offset, scale) for offset1 in offset]

        x, y = offset

        if self.do_transpose:
            if self.mirror:
                y = -y
            if self.invert:
                x = -x
        
        if scale is None:
            scale = self.scale * self.sch.node_spacing * self.size
            
        return centre + np.dot((x * self.w, y * self.h), self.R(angle_offset)) * scale


class Bipole(StretchyCpt):
    """Bipole"""

    can_mirror = True
    can_invert = True    
    can_scale = True

    node_pinnames = ('1', '2')

    pins = {'1' : ('lx', -0.5, 0),
            '2' : ('rx', 0.5, 0)}

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes[0:2]

        if self.wire:
            # With this option, draw component as a piece of wire.
            # This is useful for hiding the control voltage source
            # required for a CCVS and a CCCS.
            s = r'  \draw[-, %s] (%s) to (%s);''\n' % (
                self.args_str, n1.s, n2.s)
            return s

        tikz_cpt = self.tikz_cpt
        if self.variable:
            if self.type in ('C', 'R', 'L'):
                tikz_cpt = 'v' + tikz_cpt
            else:
                raise ValueError('Component %s not variable' % self.name)

        cpts = {'C' : {'name': 'capacitor',
                       'kinds': {'electrolytic': 'eC', 'polar': 'pC', 'variable': 'vC'},
                       'styles' : {}},
                'L' : {'name': 'inductor', 'kinds': {'variable': 'vL'},
                       'styles' : {}},
                'BAT' : {'name': 'battery', 'kinds': {'cell1': 'battery1'},
                         'styles' : {}},                
                'D' : {'name': 'diode',
                       'kinds': {'led': 'leD', 'photo': 'pD', 'schottky' : 'sD',
                                 'zener': 'zD', 'zzener': 'zzD', 'tunnel' : 'tD'},
                       'styles' : {'empty': '', 'full': '*', 'stroke': '-'}}}
        if self.type in cpts:
            cpt = cpts[self.type]
            
            if self.kind is not None:
                kinds = cpt['kinds']
                if self.kind not in kinds:
                    raise ValueError('Unknown %s kind %s: known kinds %s'
                                     % (cpt['name'], self.kind, ', '.join(kinds.keys())))
                tikz_cpt = kinds[self.kind]

            if self.style is not None:
                styles = cpt['styles']                
                if self.style not in styles:
                    raise ValueError('Unknown %s style %s: known styles %s'
                                     % (cpt['name'], self.style, ', '.join(styles.keys())))
                tikz_cpt += styles[self.style]

        elif self.type == 'MISC':
            if self.kind is None:
                raise ValueError('kind must be specified for %s' % self)
            tikz_cpt = self.kind

        label_pos = '_'
        voltage_pos = '^'
        if ((self.type in ('V', 'I', 'E', 'F', 'G', 'H', 'BAT')
             and self.sch.circuitikz_date < '2016/01/01')
            or (self.type in ('I', 'F', 'G')
                and self.sch.circuitikz_date >= '2017/05/28')):

            # Old versions of circuitikz expect the positive node
            # first, except for voltage and current sources!  So
            # swap the nodes otherwise they are drawn the wrong
            # way around.
            n1, n2 = n2, n1
                
            if self.right or self.up:
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                label_pos = '^'
                voltage_pos = '_'
        else:
            if self.left or self.down:
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                label_pos = '^'
                voltage_pos = '_'

        # Add modifier to place voltage label on other side
        # from component identifier label.
        if 'v' in self.opts:
            self.opts['v' + voltage_pos] = self.opts.pop('v')

        # Reversed voltage.
        if 'vr' in self.opts:
            self.opts['v' + voltage_pos + '>'] = self.opts.pop('vr')

        current_pos = label_pos
        # Add modifier to place current label on other side
        # from voltage marks.
        if 'i' in self.opts:
            self.opts['i' + current_pos] = self.opts.pop('i')

        # Reversed current.
        if 'ir' in self.opts:
            self.opts['i' + current_pos + '<'] = self.opts.pop('ir')

        if 'l' in self.opts:
            self.opts['l' + label_pos] = self.opts.pop('l')

        if self.type == 'O':
            node_pair_str = ''
        else:
            node_pair_str = self._node_pair_str(n1, n2, **kwargs)

        args_str = self.args_str
        args_str2 = ','.join([self.voltage_str, self.current_str, self.flow_str])

        if self.mirror:
            args_str2 += ', mirror'
        if self.invert:
            args_str2 += ', invert'            

        if self.scale != 1.0:
            args_str2 += ', bipoles/length=%.2fcm' % (self.sch.cpt_size * self.scale)

        label_str = self.label_make(label_pos, **kwargs)
            
        s = r'  \draw[%s] (%s) to [%s,%s,%s,%s,n=%s] (%s);''\n' % (
            args_str, n1.s, tikz_cpt, label_str, args_str2,
            node_pair_str, self.s, n2.s)
        return s

    
class Shape(FixedCpt):
    """General purpose shape"""

    default_aspect = 1.0
    can_mirror = True    
    can_invert = True
    pinlabels = {}

    auxiliary = {'mid' : ('c', 0.0, 0.0),
                 'bl' : ('l', -0.5, -0.5),
                 'br' : ('r', 0.5, -0.5),
                 'top' : ('t', 0, 0.5),
                 'tl' : ('l', -0.5, 0.5),
                 'tr' : ('r', 0.5, 0.5)}

    @property
    def width(self):
        return self.w * self.size * self.sch.node_spacing

    @property
    def height(self):
        return self.h * self.size * self.sch.node_spacing

    def pinpos_rotate(self, pinpos, angle):
        """Rotate pinpos by multiple of 90 degrees.  pinpos is either 'l',
        't', 'r', 'b'."""

        pin_positions = ['l', 't', 'r', 'b']
        if pinpos not in pin_positions:
            return pinpos
            
        index = pin_positions.index(pinpos)
        angle = int(angle)
        if angle < 0:
            angle += 360

        angles = (0, 90, 180, 270)
        if angle not in angles:
            raise ValueError('Cannot rotate pinpos %s by %s' % (pinpos, angle))

        index += angles.index(angle)
        pinpos = pin_positions[index % len(pin_positions)]
        return pinpos

    def parse_pinlabels(self):

        # pinlabels, pinlabels=, pinlabels=auto  label connected pins with defined labels
        # pinlabels={pin1, pin2, ...} label specified pins by pinname
        # pinlabels=all  label all pinlabels (pins connected or not)
        # pinlabels=none label no pinlabels 

        def pinlabel(nodename):

            fields = nodename.split('.')
            pinname = fields[-1]
            
            try:
                return self.pinlabels[pinname]
            except:
                pass
            return ''

        prefix = self.name + '.'

        # For backwards compatibility, check pins option.
        pinlabels = self.opts.get('pins', self.default_pins)
        if pinlabels is None:
            pinlabels = self.opts.get('pinlabels', 'none')
            
        if pinlabels == 'none':
            return {}
        elif pinlabels in ('', 'auto', 'connected'):
            return {name:pinlabel(name) for name in self.ref_node_names if pinlabel(name) != ''}
        elif pinlabels == 'all':
            return {name:pinlabel(name) for name in self.pin_node_names if pinlabel(name) != ''}
        else:
            if pinlabels[0] != '{':
                raise ValueError('Expecting { for pinlabels in %s' % self)
            if pinlabels[-1] != '}':
                raise ValueError('Expecting } for pinlabels in %s' % self)
            pinlabels = pinlabels[1:-1]
            foo = {}
            for pindef in pinlabels.split(','):
                pindef = pindef.strip()
                fields = pindef.split('=')
                if len(fields) > 1:
                    foo[prefix + fields[0].strip()] = fields[1].strip()
                else:
                    pinname = pindef
                    if pindef in self.pinlabels:
                        pinname = self.pinlabels[pindef]
                    foo[prefix + pindef] = pinname
            return foo

    def parse_pinnodes(self):

        # pinnodes, pinnodes=, pinnodes=auto  show connected pinnodes
        # pinnodes={pin1, pin2, ...} show specified pinnodes by pinname
        # pinnodes=all  show all pinnodes (connected or not)
        # pinnodes=none show no pinnodes       

        # For backwards compatibility, check anchors option.        
        pinnodes = self.opts.get('anchors', None)
        if pinnodes is None:
            pinnodes = self.opts.get('pinnodes', 'none')
            
        if pinnodes == 'none':
            return []
        elif pinnodes in ('', 'connected', 'auto'):
            return [name for name in self.ref_node_names]
        elif pinnodes == 'all':
            # Perhaps show drawing pins as well?
            return self.pin_node_names
        else:
            if pinnodes[0] != '{':
                raise ValueError('Expecting { for pinnodes in %s' % self)
            if pinnodes[-1] != '}':
                raise ValueError('Expecting } for pinnodes in %s' % self)
            pinnodes = pinnodes[1:-1]
            return [self.name + '.' + pinnode for pinnode in pinnodes.split(',')]

    def parse_pinnames(self):

        # pinnames, pinnames=, pinnames=auto  show connected pinnames
        # pinnames={pin1, pin2, ...} show specified pinnames by pinname
        # pinnames=all  show all pinnames (connected or not)
        # pinnames=none show no pinnames       

        pinnames = self.opts.get('pinnames', 'none')

        if pinnames == 'none':
            return []
        elif pinnames in ('', 'connected', 'auto'):
            return [name for name in self.ref_node_names]
        elif pinnames == 'all':
            # Perhaps show drawing pins as well?
            return self.pin_node_names
        else:
            if pinnames[0] != '{':
                raise ValueError('Expecting { for pinnames in %s' % self)
            if pinnames[-1] != '}':
                raise ValueError('Expecting } for pinnames in %s' % self)
            pinnames = pinnames[1:-1]
            return [self.name + '.' + pinname for pinname in pinnames.split(',')]

    def parse_pindefs(self):

        # pindefs={pin1=alias1, pin2=alias2, ...} define pins

        pindefs = self.opts.get('pindefs', None)
        if pindefs is None:
            return {}

        if pindefs[0] != '{':
            raise ValueError('Expecting { for pindefs in %s' % self)
        if pindefs[-1] != '}':
            raise ValueError('Expecting } for pindefs in %s' % self)
        pindefs = pindefs[1:-1]
        prefix = self.name + '.'
        foo = {}
        for pindef in pindefs.split(','):
            fields = pindef.split('=')
            if len(fields) < 2:
                raise ValueError('Expecting = in pindef %s' % pindef)
            pinname = fields[1]
            pindef = fields[0]
            if pinname not in self.pins:
                raise ValueError('Unknown pin %s in pindef' % pinname)
            foo[pinname] = pindef
        return foo
    
    def process_pinlabels(self):

        pinlabels = self.parse_pinlabels()

        for nodename, pinlabel in pinlabels.items():
            # Add pin to nodes so that it will get allocated a coord.
            node = self.sch._node_add(nodename, self, auxiliary=True)
            node.ref = self
            node.pin = not self.fakepin(node.basename)                        
            node.pinpos = self.pinpos(node.basename)

            # TODO, perhaps use pinlabel to indicate clock?
            node.clock = pinlabel != '' and pinlabel[0] == '>'
            if node.clock:
                # Remove clock designator
                pinlabel = pinlabel[1:]

            node.pinlabel = pinlabel

    def process_pinnodes(self):

        pinnodes = self.parse_pinnodes()
        for pinnode in pinnodes:
            # Add pin to nodes so that it will get allocated a coord.
            node = self.sch._node_add(pinnode, self, auxiliary=True)
            node.ref = self
            node.pin = not self.fakepin(node.basename)            
            node.pinpos = self.pinpos(node.basename)
            self.drawn_pins.append(node)

    def process_pinnames(self):

        pinnames = self.parse_pinnames()
        for pinname in pinnames:
            # Add pin to nodes so that it will get allocated a coord.
            node = self.sch._node_add(pinname, self, auxiliary=True)
            node.ref = self
            node.pin = not self.fakepin(node.basename)
            node.pinpos = self.pinpos(node.basename)
            node.pinname = node.basename
            
    def setup(self):
        super(Shape, self).setup()

        self.process_pinnodes()
        self.process_pinlabels()
        self.process_pinnames()

        # Ensure all the shape nodes are marked as pins.
        for node in self.nodes:
            node.pin = not self.fakepin(node.basename)
            
    def draw(self, **kwargs):

        if not self.check():
            return ''

        label = self.label(**kwargs)
        text_width = self.width * 0.8

        if 'image' in self.opts:
            # Override label with image
            image_filename = self.opts['image']
            ext = image_filename.split('.')[-1]
            if ext in ('tex', 'schtex', 'pgf'):

                label = r'\resizebox{%.2fcm}{!}{\input{%s}}' % (self.width,
                                                                image_filename)
            else:
                label = r'\includegraphics[width=%.2fcm]{%s}' % (self.width,
                                                                 image_filename)
            # This affects the image positioning.
            text_width = self.width

        args_str = self.args_str
        if not self.nodraw:
            args_str += ', draw'
        
        # shape border rotate rotates the box but not the text
        s = r'  \draw (%s) node[%s, thick, inner sep=0pt, minimum width=%.2fcm, minimum height=%.2fcm, text width=%.2fcm, align=center, shape border rotate=%s, %s] (%s) {%s};''\n'% (
            self.centre, self.shape, self.width, self.height, 
            text_width, self.angle, args_str, self.s, label)
        return s


class Cable(Shape):
    """Cable"""
 
    default_aspect = 4
    a = 0.3
    pins = {'in+' : ('l', -0.5, a),
            'in' : ('l', -0.5, 0),
            'in-' : ('l', -0.5, -a),
            'ignd' : ('l', -0.5, -0.45),
            'out+' : ('r', 0.5, a),
            'out' : ('r', 0.5, 0),
            'out-' : ('r', 0.5, -a),
            'ognd' : ('r', 0.475, -0.45),
            't' : ('c', 0, 0.5),                        
            'b' : ('c', 0, -0.5)}
    
    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.centre
        # Tikz uses height for length
        length = self.width
        width = self.height

        kind = self.kind
        if kind is None:
            kind = 'coax'

        if kind not in ('coax', 'twinax', 'twistedpair', 'shieldedtwistedpair',
                        'tline'):
            raise ValueError('Unknown cable kind %s' % kind)
            
        s = ''

        if kind in ('coax', 'twinax', 'shieldedtwistedpair', 'tline'):
            xscale = -1.025

            q = self.tf(centre, ((0.0125, 0)))
        
            s += r'  \draw[%s] (%s) node[cylinder, draw, rotate=%s, minimum width=%scm, minimum height=%scm, xscale=%s] {};''\n' % (self.args_str, q, self.angle, width, length, xscale)        

        if kind == 'tline':
            s += self.draw_label(centre, **kwargs)
        else:
            q = self.tf(centre, ((0, 0.9)))            
            s += self.draw_label(q, **kwargs)            
            
        if self.kind in ('twistedpair', 'shieldedtwistedpair'):
            # Needs to be even...
            twists = 4
            sections = twists * 4
            deltax = length / sections
            startx = centre.x - length / 2
            w = self.a * width
            
            x = startx
            y = self.node('mid').pos.y
            s += r'  \draw (%.2f,%.2f)' % (x, y + w)            
            # 0 -1 0 1
            for m in range(sections):
                x += deltax
                a = [0, -1, 0, 1][m % 4] * w
                s += r' %s (%.2f,%.2f)' % (('cos', 'sin')[m % 2], x, y + a)
            s += ';\n'

            x = startx
            s += r'  \draw (%.2f,%.2f)' % (x, y - w)
            # 0 1 0 -1            
            for m in range(sections):
                x += deltax
                a = [0, -1, 0, 1][m % 4] * w
                s += r' %s (%.2f,%.2f)' % (('cos', 'sin')[m % 2], x, y - a)
            s += ';\n'

        elif self.kind == 'coax':
            s += r'  \draw[-] (%s) to (%s);''\n' % (self.node('in').s,
                                                    self.node('out').s)
        elif self.kind == 'twinax':
            s += r'  \draw[-] (%s) to (%s);''\n' % (self.node('in-').s,
                                                    self.node('out-').s)
            s += r'  \draw[-] (%s) to (%s);''\n' % (self.node('in+').s,
                                                    self.node('out+').s)
        return s

    
class Transistor(FixedCpt):
    """Transistor"""
    
    can_mirror = True
    can_scale = True
    can_invert = True    

    @property
    def pins(self):
        xpins = [[self.npins, self.inpins], [self.ppins, self.ippins]]
        if self.classname in ('Qpnp', 'Mpmos', 'Jpjf'):
            return xpins[not self.mirror][self.invert]
        else:
            return xpins[self.mirror][self.invert]

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.nodes
        centre = (n1.pos + n3.pos) * 0.5

        xscale = self.scale * 2
        yscale = self.scale * 2
        if self.mirror:
            yscale = -yscale
        if self.invert:
            xscale = -xscale            
        
        s = r'  \draw (%s) node[%s, %s, xscale=%s, yscale=%s, rotate=%d] (%s) {};''\n' % (
            centre, self.tikz_cpt, self.args_str, xscale, yscale,
            self.angle, self.s)
        s += self.draw_label(centre, **kwargs)

        # Add additional wires.  These help to compensate for the
        # slight differences in sizes of the different transistors.
        if self.tikz_cpt in ('pnp', 'npn'):
            s += r'  \draw (%s.C) -- (%s) (%s.B) -- (%s) (%s.E) -- (%s);''\n' % (
                self.s, n1.s, self.s, n2.s, self.s, n3.s)
        else:
            s += r'  \draw (%s.D) -- (%s) (%s.G) -- (%s) (%s.S) -- (%s);''\n' % (
                self.s, n1.s, self.s, n2.s, self.s, n3.s)

        return s


class BJT(Transistor):
    """BJT"""

    node_pinnames = ('e', 'b', 'c')    
    ppins = {'e' : ('lx', 1, 0),
             'b' : ('lx', 0, 0.8),
             'c' : ('lx', 1, 1.6)}
    npins = {'e' : ('lx', 1, 1.6),
             'b' : ('lx', 0, 0.8),
             'c' : ('lx', 1, 0)}
    ippins = {'e' : ('lx', 0, 0),
             'b' : ('lx', 1, 0.8),
             'c' : ('lx', 0, 1.6)}
    inpins = {'e' : ('lx', 0, 1.6),
             'b' : ('lx', 1, 0.8),
             'c' : ('lx', 0, 0)}        

    
class JFET(Transistor):
    """JFET"""

    node_pinnames = ('d', 'g', 's')    
    ppins = {'d' : ('lx', 1, 0),
             'g' : ('lx', 0, 1.04),
             's' : ('lx', 1, 1.5)}
    npins = {'d' : ('lx', 1, 1.5),
             'g' : ('lx', 0, 0.46),
             's' : ('lx', 1, 0)}
    ippins = {'d' : ('lx', 0, 0),
             'g' : ('lx', 1, 1.04),
             's' : ('lx', 0, 1.5)}
    inpins = {'d' : ('lx', 0, 1.5),
             'g' : ('lx', 1, 0.46),
             's' : ('lx', 0, 0)}            


class MOSFET(Transistor):
    """MOSFET"""

    node_pinnames = ('d', 'g', 's')    
    ppins = {'d' : ('lx', 0.85, 0),
             'g' : ('lx', -0.25, 0.82),
             's' : ('lx', 0.85, 1.64)}
    npins = {'d' : ('lx', 0.85, 1.64),
             'g' : ('lx', -0.25, 0.82),
             's' : ('lx', 0.85, 0)}
    ippins = {'d' : ('lx', -0.25, 0),
             'g' : ('lx', 0.85, 0.82),
             's' : ('lx', -0.25, 1.64)}
    inpins = {'d' : ('lx', -0.25, 1.64),
             'g' : ('lx', 0.85, 0.82),
             's' : ('lx', -0.25, 0)}            


class MT(Bipole):
    """Motor"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes
        centre = (n1.pos + n2.pos) * 0.5

        s = r'  \draw (%s) node[elmech, %s, rotate=%d] (%s) {};''\n' % (
            centre, self.args_str, self.angle + 90, self.s)
        # Draw label separately, shape border rotate does not seem to work
        s += self.draw_label(centre, **kwargs)        
        s += r'  \draw (%s) |- (%s.north);''\n' % (n1.s, self.s)
        s += r'  \draw (%s.south) |- (%s);''\n' % (self.s, n2.s)
        return s

    
class MX(FixedCpt):
    """Mixer"""

    # Dubious
    can_scale = True

    node_pinnames = ('in1', 'in2', 'out')    
    pins = {'in1' : ('lx', 0.25, 0.25),
            'in2' : ('lx', -0.25, 0.25),
            'out' : ('rx', 0, 0)}

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.nodes
        centre = (n1.pos + n2.pos) * 0.5
        q = self.tf(centre, ((0, 0.35)))

        s = r'  \draw (%s) node[mixer,xscale=%s] {};''\n' % (
            centre, self.scale * self.size)
        s += self.draw_label(q, **kwargs)
        return s

    
class SP(FixedCpt):
    """Summing point"""

    # Dubious
    can_scale = True
    can_mirror = True

    node_pinnames = ('in1', 'in2', 'out', 'in3')    
    ppins = {'in1' : ('lx', -0.25, 0),
             'in2' : ('bx', 0, -0.25),
             'out' : ('rx', 0.25, 0),
             'in3' : ('tx', 0, 0.25)}
    npins = {'in1' : ('lx', -0.25, 0),
             'in2' : ('tx', 0, 0.25),
             'out' : ('rx', 0.25, 0),
             'in3' : ('bx', 0, -0.25)}    

    @property
    def pins(self):
        return self.npins if self.mirror else self.ppins
    
    def draw(self, **kwargs):

        if not self.check():
            return ''

        n = self.nodes
        centre = (n[0].pos + n[2].pos) * 0.5
        q = self.tf(centre, ((0.3, 0.3), (-0.125, 0), (0, -0.125),
                             (0, 0.125), (0, 0.125)))
        xscale = self.scale * self.size
        yscale = self.scale * self.size       
        if self.mirror:
            yscale = -yscale

        s = r'  \draw (%s) node[mixer, xscale=%s, yscale=%s, rotate=%s] {};''\n' % (
            centre, xscale, yscale, self.angle)
        s += self.draw_label(q[0], **kwargs)
        s += r'  \draw (%s) node[] {$%s$};''\n'% (q[1], self.labels[0])

        if self.mirror:
            s += r'  \draw (%s) node[] {$%s$};''\n'% (q[4], self.labels[0])
        else:
            s += r'  \draw (%s) node[] {$%s$};''\n'% (q[2], self.labels[1])
        if len(self.labels) > 2:
            s += r'  \draw (%s) node[] {$%s$};''\n'% (q[3], self.labels[2])
        return s


class SP3(SP):
    """Summing point"""

    node_pinnames = ('in1', 'in2', 'out')    
    ppins = {'in1' : ('lx', -0.25, 0),
             'in2' : ('bx', 0, -0.25),
             'out' : ('rx', 0.25, 0)}
    npins = {'in1' : ('lx', -0.25, 0),
             'in2' : ('tx', 0, 0.25),
             'out' : ('rx', 0.25, 0)}

    @property
    def pins(self):
        return self.npins if self.mirror else self.ppins


class SPpp(SP3):
    """Summing point"""

    labels = '++'


class SPpm(SP3):
    """Summing point"""

    labels = '+-'


class SPppp(SP):
    """Summing point"""

    labels = '+++'


class SPpmm(SP):
    """Summing point"""

    labels = '+--'


class SPppm(SP):
    """Summing point"""

    labels = '++-'


class TL(StretchyCpt):
    """Transmission line"""
 
    # This is dubious.  Perhaps should stretch this
    # component in proportion to size?  Applying an xscale without a
    # corresponding scale changes the ellipse.  This should be fixed
    # in circuitikz.
    can_scale = True

    node_pinnames = ('out1', 'out2', 'in1', 'in2')

    w = 1
    pins = {'in1' : ('lx', 0, 0.5),
            'in2' : ('lx', 0, 0),
            'out1' : ('rx', w, 0.5),
            'out2' : ('rx', w, 0)}
    
    @property
    def drawn_nodes(self):
        
        if self.nowires:
            return self.nodes[0:1]
        return self.nodes

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes
        centre = (n1.pos + n3.pos) * 0.5 + Pos(1 - self.w, 0)

        if self.sch.circuitikz_version >= '1.0':
            q = self.xtf(centre, ((0.23, 0),
                                  (0.20, -0.1),
                                  (-0.35, 0),
                                  (-0.27, -0.11)))
        else:
            # TODO, fix for old versions...
            q = self.xtf(centre, ((0.32, 0),
                                  (0.27, -0.145),
                                  (-0.35, 0),
                                  (-0.35, -0.145)))

        pos = n3.pos

        # Newer versions of circuitikz have a tline component with
        # wires on each end.  Rotation creates an ellipse!
        s = r'  \draw (%s) node[tlinestub, xscale=%s, rotate=%s] {};''\n' % (
            pos, self.scale, self.angle)

        s += self.draw_label(centre, **kwargs)

        # Output wire
        s += self.draw_path((q[0], n1.s))
        # Extra input wire
        s += self.draw_path((q[2], n3.s))        
        if not self.nowires:
            # Output ground wire
            s += self.draw_path((q[1], n2.s), join='|-')
            # Input ground wire            
            s += self.draw_path((q[3], n4.s), join='|-')
        return s


class TF1(FixedCpt):
    """Transformer"""

    can_scale = True    
    w = 0.8
    default_aspect = w
    node_pinnames = ('s+', 's-', 'p+', 'p-')
    pins = {'s+' : ('rx', w, 1),
            's-' : ('rx', w, 0),
            'p+' : ('lx', 0, 1),
            'p-' : ('lx', 0, 0)}
    misc = {'pdot' : (0.1 - 0.5 * w, 0.34),
            'sdot' : (0.5 * w - 0.1, 0.34),
            'link' : (0, 0.15),                 
            'label' : (0, 0.48)}

    def draw(self, link=True, **kwargs):

        if not self.check():
            return ''

        centre = self.midpoint(self.nodes[0], self.nodes[3])
        pdot_pos = self.tf(centre, self.misc['pdot'])
        sdot_pos = self.tf(centre, self.misc['sdot'])
        label_pos = self.tf(centre, self.misc['label'])

        s = ''
        if not self.nodots:
            s += r'  \draw (%s) node[circ] {};''\n' % pdot_pos
            s += r'  \draw (%s) node[circ] {};''\n' % sdot_pos
            
        s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
            label_pos, 0.5, self.s, self.label(**kwargs))

        if link:
            # TODO: allow for rotation
            width = (sdot_pos - pdot_pos).x
            link_pos = self.tf(centre, self.misc['link'])

            s += r'  \draw[<->] ([shift=(45:%.2f)]%s) arc(45:135:%.2f);''\n' % (
                width / 2, link_pos, width / 2)

        if self.classname in ('TFcore', 'TFtapcore'):
            # Draw core
            q = self.tf(centre, ((-0.05, -0.2), (-0.05, 0.2),
                                 (0.05, -0.2), (0.05, 0.2)))
            s += self.draw_path(q[0:2], style='thick')
            s += self.draw_path(q[2:4], style='thick')

        return s


class Transformer(TF1):
    """Transformer"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes[0:4]

        s = r' \draw[%s] (%s) to [inductor, scale=%s] (%s);''\n' % (
            self.args_str, n3.s, self.scale, n4.s)
        s += r' \draw[%s] (%s) to [inductor, scale=%s] (%s);''\n' % (
            self.args_str, n2.s, self.scale, n1.s)

        s += super(Transformer, self).draw(link=False, **kwargs)
        return s


class TFtap(TF1):
    """Transformer"""

    node_pinnames = ('s+', 's-', 'p+', 'p-', 'ptap', 'stap')
    w = 0.8
    pins = {'s+' : ('rx', w, 1),
            's-' : ('rx', w, 0),
            'p+' : ('lx', 0, 1),
            'p-' : ('lx', 0, 0),
            'ptap' : ('lx', 0, 0.5),
            'stap' : ('rx', w, 0.5)}
    
    @property
    def drawn_nodes(self):
        # Do not draw the taps.
        return self.nodes[0:4]

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes[0:4]

        s = r'  \draw (%s) to [inductor] (%s);''\n' % (n3.s, n4.s)
        s += r'  \draw (%s) to [inductor] (%s);''\n' % (n2.s, n1.s)

        s += super(TFtap, self).draw(link=False, **kwargs)
        return s


class K(TF1):
    """Mutual coupling"""

    def __init__(self, sch, namespace, defname, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        self.Lname1 = args[0]
        self.Lname2 = args[1]
        super(K, self).__init__(sch, namespace, defname, name,
                                cpt_type, cpt_id, string, opts_string,
                                node_names, keyword, *args[2:])

    @property
    def coords(self):
        return ((0.5, 1), (0.5, 0), (0, 1), (0, 0))

    @property
    def scales(self):
        return [self.scale] * len(self.coords)
        
    @property
    def nodes(self):

        # This needs to be more robust; currently it depends on the order
        # that the inductors are defined.
        
        # L1 and L2 need to be previously defined so we can find their nodes.
        L1 = self.sch.elements[self.Lname1]
        L2 = self.sch.elements[self.Lname2]
        # L1 is on the left; L2 is on the right
        nodes = [self.sch.nodes[n] for n in L2.node_names + L1.node_names]
        return nodes


class Gyrator(FixedCpt):
    """Gyrator"""

    node_pinnames = ('out+', 'out-', 'in+', 'in-')    
    pins = {'out+' : ('rx', 1, 1),
            'out-' : ('rx', 1, 0),
            'in+' : ('lx', 0, 1),
            'in-' : ('lx', 0, 0)}
    
    def draw(self, **kwargs):

        if not self.check():
            return ''

        yscale = self.scale
        if not self.mirror:
            yscale = -yscale

        s = r'  \draw (%s) node[gyrator, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            self.midpoint(self.nodes[0], self.nodes[3]),
            self.args_str, 0.95 * self.scale, 0.89 * yscale,
            -self.angle, self.s)        

        s += self.draw_label(self.centre, **kwargs)
        return s

        
class Potentiometer(Bipole):
    """Potentiometer  Np, Nm, No"""

    # This is not really a bipole but circuitikz treats it as such

    can_stretch = False

    node_pinnames = ('1', '2', 'wiper')    
    pins = {'1' : ('rx', 0, 0),
            '2' : ('rx', 1, 0),
            'wiper' : ('lx', 0.5, 0.3)}
    
    
class VCS(Bipole):
    """Voltage controlled source"""

    @property
    def required_node_names(self):
        return self.node_names[0:2]


class CCS(Bipole):
    """Current controlled source"""

    @property
    def required_node_names(self):
        return self.node_names[0:2]    


class SPDT(FixedCpt):
    """SPDT switch"""

    can_mirror = True
    can_invert = True    

    node_pinnames = ('1', '2', 'common')    
    pins = {'1' : ('lx', 0, 0.169),
            '2' : ('rx', 0.632, 0.338),
            'common' : ('lx', 0.632, 0)}
    
    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3 = self.nodes
        centre = n1.pos * 0.5 + (n2.pos + n3.pos) * 0.25

        args_str = 'rotate=%d' % self.angle
        
        if self.mirror:
            args_str += ', mirror'
        if self.invert:
            args_str += ', invert'            
        
        s = r'  \draw (%s) node[spdt, %s, %s] (%s) {};''\n' % (
            centre, self.args_str, args_str, self.s)            
            
        # TODO, fix label position.
        centre = (n1.pos + n3.pos) * 0.5 + Pos(0.5, -0.5)
        s += self.draw_label(centre, **kwargs)
        return s


class Box2(Shape):
    """Square box,  A rectangle is created by defining aspect."""

    shape = 'rectangle'
    pins = {'w' : ('l', -0.5, 0),
            'e' : ('r', 0.5, 0)}
    
    
class Box4(Shape):
    """Box4"""

    shape = 'rectangle'
    pins = {'w' : ('l', -0.5, 0),
            's' : ('b', 0, -0.5),
            'e' : ('r', 0.5, 0),
            'n' : ('t', 0, 0.5)}
    

class Box12(Shape):
    """Box12"""

    shape = 'rectangle'
    pins = {'wnw' : ('l', -0.5, 0.25),
            'w' : ('l', -0.5, 0),
            'wsw' : ('l', -0.5, -0.25),
            'ssw' : ('b', -0.25, -0.5),                              
            's' : ('b', 0, -0.5),
            'sse' : ('b', 0.25, -0.5),
            'ese' : ('r', 0.5, -0.25),                              
            'e' : ('r', 0.5, 0),
            'ene' : ('r', 0.5, 0.25),
            'nne' : ('t', 0.25, 0.5),
            'n' : ('t', 0, 0.5),
            'nnw' : ('t', -0.25, 0.5)}


class Box(Shape):
    """Box"""

    shape = 'rectangle'    
    pins = {'nw' : ('t', -0.5, 0.5), 'wnw' : ('l', -0.5, 0.25),
            'w' : ('l', -0.5, 0), 'wsw' : ('l', -0.5, -0.25), 
            'sw' : ('b', -0.5, -0.5), 'ssw' : ('b', -0.25, -0.5),
            's' : ('b', 0, -0.5), 'sse' : ('b', 0.25, -0.5),
            'se' : ('b', 0.5, -0.5), 'ese' : ('r', 0.5, -0.25),
            'e' : ('r', 0.5, 0), 'ene' : ('r', 0.5, 0.25),
            'ne' : ('t', 0.5, 0.5), 'nne' : ('t', 0.25, 0.5),
            'n' : ('t', 0, 0.5), 'nnw' : ('t', -0.25, 0.5)}


class Ellipse(Shape):
    """Ellipse"""

    # Ellipse needs the tikz shapes library.
    shape = 'ellipse'
    pins = {'nw' : ('t', -0.3536, 0.3536), 'wnw' : ('l', -0.4619, 0.1913),
            'w' : ('l', -0.5, 0), 'wsw' : ('l', -0.4619, -0.1913), 
            'sw' : ('b', -0.3536, -0.3536), 'ssw' : ('b', -0.1913, -0.4619),
            's' : ('b', 0, -0.5), 'sse' : ('b', 0.1913, -0.4619),
            'se' : ('r', 0.3536, -0.3536), 'ese' : ('r', 0.4619, -0.1913),
            'e' : ('r', 0.5, 0), 'ene' : ('r', 0.4619, 0.1913),
            'ne' : ('r', 0.3536, 0.35365), 'nne' : ('t', 0.1913, 0.4619),
            'n' : ('t', 0, 0.5), 'nnw' : ('t', -0.1913, 0.4619)}


class Circle(Ellipse):
    """Circle"""

    shape = 'circle'


class Circle2(Shape):
    """Circle"""

    shape = 'circle'
    pins = {'w' : ('l', -0.5, 0),
            'e' : ('r', 0.5, 0)}
    

class Circle4(Shape):
    """Circle4"""

    shape = 'circle'
    pins = {'w' : ('l', -0.5, 0),
            's' : ('b', 0, -0.5),
            'e' : ('r', 0.5, 0),
            'n' : ('t', 0, 0.5)}


class Triangle(Shape):
    """Equilateral triangle, The triangle shape can be altered by defining
    aspect."""    

    shape = 'triangle'
    
    # 1 / sqrt(3) approx 0.5774, 1 / (2 * sqrt(3)) approx 0.2887
    pins = {'n' : ('t', 0.0, 0.5774),
            'ne' : ('r', 0.25, 0.14435),
            'nw' : ('l', -0.25, 0.14435),
            'w' : ('l', -0.5, -0.2887),
            'e' : ('r', 0.5, -0.2887),
            's' : ('b', 0.0, -0.2887),
            'se' : ('b', 0.25, -0.2887),
            'sw' : ('b', -0.25, -0.2887),
            'ssw' : ('b', -0.125, -0.2887),
            'sse' : ('b', 0.125, -0.2887),
            'nne' : ('r', 0.125, 0.355),
            'nnw' : ('l', -0.125, 0.355),
            'wsw' : ('b', -0.375, -0.2887),
            'ese' : ('b', 0.375, -0.2887),
            'ene' : ('r', 0.375, -0.075),
            'wnw' : ('l', -0.375, -0.075)}

    auxiliary = {'mid' : ('c', 0.0, 0.0),
                 'bl' : ('l', -0.5, -0.2887),
                 'br' : ('r', 0.5, -0.2887),
                 'top' : ('t', 0, 0.5774),                    
                 'tl' : ('l', -0.5, 0.5774),
                 'tr' : ('r', 0.5, 0.5774)}

    def draw(self, **kwargs):

        if not self.check():
            return ''

        s = self.draw_path([self.node('top').pos, self.node('bl').pos,
                            self.node('br').pos], closed=True, style='thick')
        s += self.draw_label(self.node('mid').pos, **kwargs)

        return s

    
class TwoPort(Shape):
    """Two-port"""

    default_width = 2
    default_aspect = 1
    shape = 'rectangle'
    pins = {'wnw' : ('l', -0.5, 0.375),
            'w' : ('l', -0.5, 0),
            'wsw' : ('l', -0.5, -0.375),
            'ssw' : ('b', -0.25, -0.5),                              
            's' : ('b', 0, -0.5),
            'sse' : ('b', 0.25, -0.5),
            'ese' : ('r', 0.5, -0.375),                              
            'e' : ('r', 0.5, 0),
            'ene' : ('r', 0.5, 0.375),
            'nne' : ('t', 0.25, 0.5),
            'n' : ('t', 0, 0.5),
            'nnw' : ('t', -0.25, 0.5)}
    node_pinnames = ('ene', 'ese', 'wnw', 'wsw')
    

class TR(Box2):
    """Transfer function"""

    default_width = 1.5
    default_aspect = 1.5
    node_pinnames = ('w', 'e')    

    
class Chip(Shape):
    """General purpose chip"""

    do_transpose = True
    default_width = 2.0

    # Could allow can_scale but not a lot of point since nodes
    # will not be on the boundary of the chip.

    # TODO, tweak coord if pin name starts with \ using pinpos to
    # accomodate inverting circle.  This will require stripping of the
    # \ from the label. Alternatively, do not use inverting circle and
    # add overline to symbol name when printing.

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5))

    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.node('mid')

        s = ''
        if isinstance(self.path, list):
            # For multiple paths, e.g., isoamp
            for path in self.path:
                q = self.tf(centre.pos, path)
                s += self.draw_path(q, closed=True, style='thick')
        else:
            q = self.tf(centre.pos, self.path)
            s = self.draw_path(q, closed=True, style='thick')

        label = self.label(**kwargs)
        if label != '':
            s += r'  \draw (%s) node[text width=%.2fcm, align=center, %s] {%s};''\n'% (
                centre.s, self.width - 0.5, self.args_str, label)

        # Draw clock symbols
        for m, n in enumerate(self.nodes):
            if n.clock:
                q = self.tf(n.pos, ((0, 0.125 * 0.707), (0.125, 0), 
                                    (0, -0.125 * 0.707)))
                s += self.draw_path(q[0:3], style='thick')
        return s


class Uchip1313(Chip):
    """Chip of size 1 3 1 3"""

    pins = {'in' : ('l', -0.5, 0),
            'b1' : ('b', -0.25, -0.5),
            'vss' : ('b', 0, -0.5),
            'b3': ('b', 0.25, -0.5),
            'out': ('r', 0.5, 0),
            't1' : ('t', -0.25, 0.5),
            'vdd' : ('t', 0, 0.5),
            't3': ('t', 0.25, 0.5)}            
    

class Uchip2121(Chip):
    """Chip of size 2 1 2 1"""

    pins = {'l1' : ('l', -0.5, 0.25),
            'l2' : ('l', -0.5, -0.25),
            'vss' : ('b', 0, -0.5),
            'r2': ('r', 0.5, -0.25),
            'r1': ('r', 0.5, 0.25),
            'vdd': ('t', 0, 0.5)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD'}


class Uchip3131(Chip):
    """Chip of size 3 1 3 1"""

    pins = {'l1' : ('l', -0.5, 0.25),
            'l2' : ('l', -0.5, 0),               
            'l3' : ('l', -0.5, -0.25),
            'vss' : ('b', 0.0, -0.5),
            'r3': ('r', 0.5, -0.25),
            'r2': ('r', 0.5, 0),               
            'r1': ('r', 0.5, 0.25),
            'vdd': ('t', 0.0, 0.5)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD'}

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5))


class Uchip3333(Chip):
    """Chip of size 3 3 3 3"""

    pins = {'l1' : ('l', -0.5, 0.25),
            'l2' : ('l', -0.5, 0),               
            'l3' : ('l', -0.5, -0.25),
            'b1' : ('b', -0.25, -0.5),            
            'vss' : ('b', 0.0, -0.5),
            'b3' : ('b', 0.25, -0.5),                        
            'r3': ('r', 0.5, -0.25),
            'r2': ('r', 0.5, 0),               
            'r1': ('r', 0.5, 0.25),
            't1' : ('t', -0.25, 0.5),                        
            'vdd': ('t', 0.0, 0.5),
            't3' : ('t', 0.25, 0.5)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD'}

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5))    


class Uchip2222(Chip):
    """Chip of size 2 2 2 2"""

    pins = {'l1' : ('l', -0.5, 0.25),
            'l2' : ('l', -0.5, -0.25),
            'b1' : ('b', -0.25, -0.5),
            'b2' : ('b', 0.25, -0.5),
            'r1' : ('r', 0.5, 0.25),
            'r2' : ('r', 0.5, -0.25),
            't1' : ('t', -0.25, 0.5),
            't2' : ('t', 0.25, 0.5)}
    

class Uchip4141(Chip):
    """Chip of size 4 1 4 1"""

    pins = {'l1' : ('l', -0.5, 0.375),
            'l2' : ('l', -0.5, 0.125),
            'l3' : ('l', -0.5, -0.125),               
            'l4' : ('l', -0.5, -0.375),
            'vss' : ('b', 0.0, -0.5),
            'r4': ('r', 0.5, -0.375),
            'r3': ('r', 0.5, -0.125),               
            'r2': ('r', 0.5, 0.125),
            'r1': ('r', 0.5, 0.375),
            'vdd': ('t', 0.0, 0.5)}

class Uchip4444(Chip):
    """Chip of size 4 4 4 4"""

    pins = {'l1' : ('l', -0.5, 0.375),
            'l2' : ('l', -0.5, 0.125),
            'l3' : ('l', -0.5, -0.125),               
            'l4' : ('l', -0.5, -0.375),
            'b1' : ('b', -0.375, -0.5),
            'b2' : ('b', -0.125, -0.5),
            'b3' : ('b', .125, -0.5),
            'b4' : ('b', 0.375, -0.5),            
            'r4': ('r', 0.5, -0.375),
            'r3': ('r', 0.5, -0.125),               
            'r2': ('r', 0.5, 0.125),
            'r1': ('r', 0.5, 0.375),
            't1' : ('t', -0.375, 0.5),
            't2' : ('t', -0.125, 0.5),
            't3' : ('t', .125, 0.5),
            't4' : ('t', 0.375, 0.5)}
            

class Uchip8181(Chip):
    """Chip of size 8 1 8 1"""

    pins = {'l1' : ('l', -0.5, 0.4375),
            'l2' : ('l', -0.5, 0.3125),
            'l3' : ('l', -0.5, 0.1875),
            'l4' : ('l', -0.5, 0.0625),
            'l5' : ('l', -0.5, -0.0625),                        
            'l6' : ('l', -0.5, -0.1875),               
            'l7' : ('l', -0.5, -0.3125),
            'l8' : ('l', -0.5, -0.4375),        
            'vss' : ('b', 0.0, -0.5),
            'r8': ('r', 0.5, -0.4375),
            'r7': ('r', 0.5, -0.3125),
            'r6': ('r', 0.5, -0.1875),
            'r5': ('r', 0.5, -0.0625),               
            'r4': ('r', 0.5, 0.0625),
            'r3': ('r', 0.5, 0.1875),
            'r2': ('r', 0.5, 0.3125),            
            'r1': ('r', 0.5, 0.4375),            
            'vdd': ('t', 0.0, 0.5)}

class Uchip8888(Chip):
    """Chip of size 8 8 8 8"""

    pins = {'l1' : ('l', -0.5, 0.4375),
            'l2' : ('l', -0.5, 0.3125),
            'l3' : ('l', -0.5, 0.1875),
            'l4' : ('l', -0.5, 0.0625),
            'l5' : ('l', -0.5, -0.0625),                        
            'l6' : ('l', -0.5, -0.1875),               
            'l7' : ('l', -0.5, -0.3125),
            'l8' : ('l', -0.5, -0.4375),
            'b1' : ('b', -0.4375, -0.5),
            'b2' : ('b', -0.3125, -0.5),
            'b3' : ('b', -0.1875, -0.5),
            'b4' : ('b', -0.0625, -0.5),
            'b5' : ('b', 0.0625, -0.5),                        
            'b6' : ('b', 0.1875, -0.5),               
            'b7' : ('b', 0.3125, -0.5),
            'b8' : ('b', -0.4375, -0.5),                    
            'r8': ('r', 0.5, -0.4375),
            'r7': ('r', 0.5, -0.3125),
            'r6': ('r', 0.5, -0.1875),
            'r5': ('r', 0.5, -0.0625),               
            'r4': ('r', 0.5, 0.0625),
            'r3': ('r', 0.5, 0.1875),
            'r2': ('r', 0.5, 0.3125),            
            'r1': ('r', 0.5, 0.4375),
            't1' : ('t', -0.4375, 0.5),
            't2' : ('t', -0.3125, 0.5),
            't3' : ('t', -0.1875, 0.5),
            't4' : ('t', -0.0625, 0.5),
            't5' : ('t', 0.0625, 0.5),                        
            't6' : ('t', 0.1875, 0.5),               
            't7' : ('t', 0.3125, 0.5),
            't8' : ('t', -0.4375, 0.5)}

class Mux(Chip):
    """Multiplexer"""

    @property
    def path(self):
        return ((-0.25, -0.5), (-0.25, 0.5), (0.25, 0.25), (0.25, -0.25), (-0.25, -0.5))    

    
class Umux21(Mux):
    """Multiplexer 2 to 1"""

    pins = {'l1' : ('l', -0.25, 0.25),
            'l2' : ('l', -0.25, -0.25),
            'b' : ('b', 0, -0.375),
            't' : ('t', 0, 0.375),            
            'r': ('r', 0.25, 0.0)}


class Umux41(Mux):
    """Multiplexer 4 to 1"""

    pins = {'l1' : ('l', -0.25, 0.375),
            'l2' : ('l', -0.25, 0.125),
            'l3' : ('l', -0.25, -0.125),
            'l4' : ('l', -0.25, -0.375),            
            'b1' : ('b', -0.125, -0.4375),
            'b2' : ('b', 0.125, -0.3125),
            't1' : ('t', -0.125, 0.4375),
            't2' : ('t', 0.125, 0.3125),                        
            'r': ('r', 0.25, 0)}

    
class Umux42(Mux):
    """Multiplexer 4 to 2"""

    pins = {'l1' : ('l', -0.25, 0.375),
            'l2' : ('l', -0.25, 0.125),
            'l3' : ('l', -0.25, -0.125),
            'l4' : ('l', -0.25, -0.375),            
            'b1' : ('b', -0.125, -0.4375),
            'b2' : ('b', 0.125, -0.3125),
            't1' : ('t', -0.125, 0.4375),
            't2' : ('t', 0.125, 0.3125),                        
            'r1': ('r', 0.25, 0.125),
            'r2': ('r', 0.25, -0.125)}
    
    
class Uadc(Chip):
    """ADC"""

    pins = {'in' : ('l', -0.5, 0),
            'in+' : ('l', -0.4375, 0.125),
            'in-' : ('l', -0.4375, -0.125),
            'vref-' : ('l', -0.375, -0.25),
            'vref+' : ('l', -0.375, 0.25),                              
            'avss' : ('b', -0.1, -0.5),
            'vss' : ('b', 0.1, -0.5),            
            'dvss' : ('b', 0.3, -0.5),
            'clk' : ('r', 0.5, -0.25),
            'data' : ('r', 0.5, 0),
            'fs' : ('r', 0.5, 0.25),
            'dvdd' : ('t', 0.3, 0.5),
            'vdd' : ('t', 0.1, 0.5),            
            'avdd' : ('t', -0.1, 0.5)}
    
    pinlabels = {'vref-' : 'VREF-', 'vref+' : 'VREF+',
                 'vss': 'VSS', 'vdd' : 'VDD',
                 'dvss': 'DVSS', 'dvdd' : 'DVDD',                 
                 'avss': 'AVSS', 'avdd' : 'AVDD',                 
                 'clk' : '>', 'data' : 'DATA', 'fs' : 'FS'}

    @property
    def path(self):
        return ((-0.5, 0.0), (-0.25, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.25, 0.5))

    
class Uregulator(Chip):
    """Voltage regulator"""

    default_aspect = 4.0 / 3.0
    pins = {'in' : ('l', -0.5, 0),
            'en' : ('b', -0.25, -0.5),
            'gnd' : ('b', 0, -0.5),               
            'out': ('r', 0.5, 0)}
    
    pinlabels = {'en' : 'E', 'gnd' : 'GND'}
    

class Udac(Chip):
    """DAC"""

    pins = {'out' : ('r', 0.5, 0),
            'out+' : ('r', 0.4375, 0.125),
            'out-' : ('r', 0.4375, -0.125),               
            'vref-' : ('r', 0.375, -0.25),
            'vref+' : ('r', 0.375, 0.25),
            'avss' : ('b', 0.1, -0.5),
            'vss' : ('b', -0.1, -0.5),
            'dvss' : ('b', -0.3, -0.5),
            'clk' : ('l', -0.5, -0.25),
            'data' : ('l', -0.5, 0),
            'fs' : ('l', -0.5, 0.25),
            'dvdd' : ('t', -0.3, 0.5),
            'vdd' : ('t', -0.1, 0.5),            
            'avdd' : ('t', 0.1, 0.5)}
    
    pinlabels = {'vref-' : 'VREF-', 'vref+' : 'VREF+',
                 'vss': 'VSS', 'vdd' : 'VDD',
                 'dvss': 'DVSS', 'dvdd' : 'DVDD',                 
                 'avss': 'AVSS', 'avdd' : 'AVDD',
                 'clk' : '>', 'data' : 'DATA', 'fs' : 'FS'}
    
    @property
    def path(self):
        return ((-0.5, -0.5), (0.25, -0.5), (0.5, 0), (0.25, 0.5), (-0.5, 0.5))


class Udiffamp(Chip):
    """Differential amplifier.  It is not automatically annotated with + and - symbols for inputs.
    This can be achieved using pinlabels={in+, in-}."""

    default_width = 1.0

    pins = {'in+' : ('l', -0.5, 0.25),
            'in-' : ('l', -0.5, -0.25),
            'vss' : ('b', 0, -0.25),
            'out' : ('r', 0.5, 0),
            'vdd' : ('t', 0, 0.25)}

    pinlabels = {'in+' : '$+$', 'in-' : '$-$', 'vss' : 'VSS', 'vdd' : 'VDD'}

    @property
    def path(self):
        return ((-0.5, 0.5), (-0.5, -0.5), (0.5, 0))


class Ubuffer(Chip):
    """Buffer with power supplies"""

    default_width = 1.0

    pins = {'in' : ('l', -0.5, 0),
            'vss' : ('b', 0, -0.25),
            'out' : ('r', 0.5, 0),
            'vdd': ('t', 0, 0.25),
            'vdd1': ('t', -0.25, 0.375),
            'vdd2': ('t', 0.25, 0.125),                        
            'en': ('b', -0.25, -0.375)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD', 'en' : 'E'}

    @property
    def path(self):
        return ((-0.5, 0.5), (0.5, 0), (-0.5, -0.5))


class Uinverter(Chip):
    """Inverter with power supplies"""

    default_width = 1.0

    pins = {'in' : ('l', -0.5, 0),
            'vss' : ('b', 0, -0.22),
            'out' : ('r', 0.5, 0),
            'vdd': ('t', 0, 0.22),
            'en': ('b', -0.25, -0.37)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD', 'en' : 'E'}

    @property
    def path(self):
        w = 0.05
        return ((-0.5, 0.5), (0.5 - 2 * w, 0), (-0.5, -0.5))

    def draw(self, **kwargs):

        s = super(Uinverter, self).draw(**kwargs)

        # Append inverting circle.
        centre = self.node('mid')                
        q = self.tf(centre.pos, ((0.45, 0)))
        s += r'  \draw[thick] (%s) node[ocirc, scale=%s, %s] {};''\n' % (
            q, 1.8 * self.size * self.scale, self.args_str)
        return s


class Udiffdriver(Chip):
    """Differential driver with power supplies"""

    default_width = 1.0

    pins = {'in' : ('l', -0.5, 0),
            'vss' : ('b', -0.15, -0.3),
            'out+' : ('r', -0.02, 0.25),
            'out-' : ('r', 0.1, -0.25),            
            'vdd': ('t', 0, 0.22),
            'en': ('b', -0.3, -0.4)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD', 'en' : 'E'}

    @property
    def path(self):
        w = 0.05
        return ((-0.5, 0.5), (0.5 - 2 * w, 0), (-0.5, -0.5))

    def draw(self, **kwargs):

        s = super(Udiffdriver, self).draw(**kwargs)

        # Append inverting circle.
        centre = self.node('mid')                
        q = self.tf(centre.pos, ((0.05, -0.25)))
        s += r'  \draw[thick] (%s) node[ocirc, scale=%s, %s] {};''\n' % (
            q, 1.8 * self.size * self.scale, self.args_str)
        return s    

    
class Flipflop(Chip):

    default_width = 1.0    

    
class Udff(Flipflop):
    """D flip-flop"""

    default_pins = '{d, clk, q}'
    
    pins = {'d' : ('l', -0.5, 0.25),
            'clk' : ('l', -0.5, 0),               
            'vss' : ('b', 0.0, -0.5),
            '/q': ('r', 0.5, -0.25),
            'q': ('r', 0.5, 0.25),
            'vdd': ('t', 0.0, 0.5)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD',
                 'd' : 'D', 'q' : 'Q', '/q' : '$\overline{\mathrm{Q}}$', 'clk' : '>'}
    

class Ujkff(Flipflop):
    """JK flip-flop"""

    default_pins = '{j, k, clk, q}'
    
    pins = {'j' : ('l', -0.5, 0.25),
            'clk' : ('l', -0.5, 0),
            'k' : ('l', -0.5, -0.25),               
            'vss' : ('b', 0.0, -0.5),
            '/q': ('r', 0.5, -0.25),
            'q': ('r', 0.5, 0.25),
            'vdd': ('t', 0.0, 0.5)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD',
                 'j' : 'J', 'k' : 'K',
                 'q' : 'Q', '/q' : '$\overline{\mathrm{Q}}$', 'clk' : '>'}    


class Urslatch(Flipflop):
    """RS latch"""

    default_pins = '{r, s, q}'
    
    pins = {'r' : ('l', -0.5, 0.25),
            's' : ('l', -0.5, -0.25),               
            'vss' : ('b', 0.0, -0.5),
            '/q': ('r', 0.5, -0.25),
            'q': ('r', 0.5, 0.25),
            'vdd': ('t', 0.0, 0.5)}
    
    pinlabels = {'vss' : 'VSS', 'vdd' : 'VDD',
                 'r' : 'R', 's' : 'S',
                 'q' : 'Q', '/q' : '$\overline{\mathrm{Q}}$'}    
    
    
class Eopamp(Chip):
    """This is for an opamp created with the E netlist type as used for
    simulation."""

    can_scale = True
    can_mirror = True
    do_transpose = False
    default_width = 1.0

    # The Nm node is not used (ground).
    node_pinnames = ('out', '', 'in+', 'in-')
    
    ppins = {'out' : ('rx', 1.25, 0.0),
             'in+' : ('lx', -1.25, 0.5),
             'in-' : ('lx', -1.25, -0.5),
             'vdd' : ('t', 0, 0.5),
             'vdd2' : ('t', -0.45, 0.755),
             'vss2' : ('b', -0.45, -0.755),
             'vss' : ('b', 0, -0.5),
             'ref' : ('b', 0.45, -0.245),
             'r+' : ('l', -0.85, 0.25),
             'r-' : ('l', -0.85, -0.25)}

    npins = {'out' : ('rx', 1.25, 0.0),
             'in-' : ('lx', -1.25, 0.5),
             'in+' : ('lx', -1.25, -0.5),
             'vdd' : ('t', 0, 0.5),
             'vdd2' : ('t', -0.45, 0.755),
             'vss2' : ('b', -0.45, -0.755),
             'vss' : ('b', 0, -0.5),
             'ref' : ('b', 0.45, -0.245),
             'r-' : ('l', -0.85, 0.25),
             'r+' : ('l', -0.85, -0.25)}    

    pinlabels = {'vdd' : 'VDD', 'vss' : 'VSS'}

    @property
    def pins(self):
        return self.npins if (self.mirrorinputs ^ self.mirror) else self.ppins

    def draw(self, **kwargs):

        if not self.check():
            return ''

        yscale = 2 * 0.95 * self.scale        
        if not (self.mirror ^ self.mirrorinputs):
            yscale = -yscale

        centre = self.node('mid')

        # The circuitikz opamp has short wires on input and output.
        
        # Note, scale scales by area, xscale and yscale scale by length.
        s = r'  \draw (%s) node[op amp, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre.s,
            self.args_str, 2 * self.scale * 0.95, yscale,
            -self.angle, self.s)
        if not self.nowires:
            s += r'  \draw (%s.out) |- (%s);''\n' % (self.s, self.node('out').s)
            s += r'  \draw (%s.+) |- (%s);''\n' % (self.s, self.node('in+').s)
            s += r'  \draw (%s.-) |- (%s);''\n' % (self.s, self.node('in-').s)
        s += self.draw_label(centre.s, **kwargs)
        return s


class Efdopamp(Chip):
    """This is for a fully differential opamp created with the E netlist
    type.  See also Ufdopamp for a fully differential opamp created
    with the U netlist type."""    

    can_scale = True
    can_mirror = True
    do_transpose = False
    default_width = 1.0    

    node_pinnames = ('out+', 'out-', 'in+', 'in-')

    ppins = {'out+' : ('r', 0.85, -0.5),
             'out-' : ('r', 0.85, 0.5),                
             'in+' : ('l', -1.25, 0.5),
             'vocm' : ('l', -0.85, 0),                                       
             'in-' : ('l', -1.25, -0.5),
             'vdd' : ('t', -0.25, 0.645),
             'vss' : ('b', -0.25, -0.645),
             'r+' : ('l', -0.85, 0.25),
             'r-' : ('l', -0.85, -0.25)}

    npins = {'out-' : ('r', 0.85, -0.5),
             'out+' : ('r', 0.85, 0.5),                
             'in-' : ('l', -1.25, 0.5),
             'vocm' : ('l', -0.85, 0),                          
             'in+' : ('l', -1.25, -0.5),
             'vdd' : ('t', -0.25, 0.645),
             'vss' : ('b', -0.25, -0.645),
             'r-' : ('l', -0.85, 0.25),
             'r+' : ('l', -0.85, -0.25)}    
    
    pinlabels = {'vdd' : 'VDD', 'vss' : 'VSS'}
    
    @property
    def pins(self):
        return self.npins if (self.mirrorinputs ^ self.mirror) else self.ppins
    
    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.node('mid')

        yscale = 2 * 0.952 * self.scale
        if not (self.mirror ^ self.mirrorinputs):        
            yscale = -yscale

        s = r'  \draw (%s) node[fd op amp, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre.s, self.args_str, 2 * self.scale * 0.95, yscale,
            -self.angle, self.s)
        s += r'  \draw (%s.out +) |- (%s);''\n' % (self.s, self.node('out+').s)
        s += r'  \draw (%s.out -) |- (%s);''\n' % (self.s, self.node('out-').s)
        s += r'  \draw (%s.+) |- (%s);''\n' % (self.s, self.node('in+').s)
        s += r'  \draw (%s.-) |- (%s);''\n' % (self.s, self.node('in-').s)  

        s += self.draw_label(centre.s, **kwargs)
        return s

class Eamp(Chip):
    """Amplifier."""

    can_scale = True
    do_transpose = False    
    default_width = 1.0

    # The Nm and Ncm nodes are not used (ground).
    node_pinnames = ('out', '', 'in', '')
    
    pins = {'out' : ('rx', 1.25, 0.0),
            'in' : ('lx', -1.25, 0.0),
            'vdd' : ('t', 0, 0.5),
            'vdd2' : ('t', -0.45, 0.755),
            'vss2' : ('b', -0.45, -0.755),
            'vss' : ('b', 0, -0.5),
            'ref' : ('b', 0.45, -0.245),
            'r+' : ('l', -0.85, 0.25),
            'r-' : ('l', -0.85, -0.25)}

    pinlabels = {'vdd' : 'VDD', 'vss' : 'VSS'}


    def draw(self, **kwargs):

        if not self.check():
            return ''

        centre = self.node('mid')
        yscale = 2 * 0.952 * self.scale        

        # The circuitikz buffer has short wires on input and output.
        
        # Note, scale scales by area, xscale and yscale scale by length.
        s = r'  \draw (%s) node[buffer, %s, xscale=%.3f, yscale=%.3f, rotate=%d] (%s) {};''\n' % (
            centre.s,
            self.args_str, 2 * self.scale * 0.95, yscale,
            -self.angle, self.s)
        if not self.nowires:
            s += r'  \draw (%s.out) |- (%s);''\n' % (self.s, self.node('out').s)
            s += r'  \draw (%s.-) |- (%s);''\n' % (self.s, self.node('in').s)
        s += self.draw_label(centre.s, **kwargs)
        return s


class Uopamp(Chip):
    """This is for an opamp created with the U netlist type.  It has no wires.
    See also Opamp for an opamp created with the E netlist type."""

    can_mirror = True

    auxiliary = {'lin+' : ('c', -0.375, 0.25),                        
                 'lin-' : ('c', -0.375, -0.25)}
    auxiliary.update(Chip.auxiliary)
    
    ppins = {'out' : ('r', 0.5, 0.0),
             'in+' : ('l', -0.5, 0.25),
             'in-' : ('l', -0.5, -0.25),
             'vdd' : ('t', 0, 0.25),
             'vdd2' : ('t', -0.225, 0.365),
             'vss2' : ('b', -0.225, -0.365),
             'vss' : ('b', 0, -0.25),
             'ref' : ('b', 0.225, -0.135),
             'r+' : ('l', -0.5, 0.125),
             'r-' : ('l', -0.5, -0.125)}

    npins = {'out' : ('r', 0.5, 0.0),
             'in-' : ('l', -0.5, 0.25),
             'in+' : ('l', -0.5, -0.25),
             'vdd' : ('t', 0, 0.25),
             'vdd2' : ('t', -0.225, 0.365),
             'vss2' : ('b', -0.225, -0.365),
             'vss' : ('b', 0, -0.25),
             'ref' : ('b', 0.225, -0.135),
             'r-' : ('l', -0.5, 0.125),
             'r+' : ('l', -0.5, -0.125)}    

    pinlabels = {'vdd' : 'VDD', 'vss' : 'VSS'}

    @property
    def pins(self):
        return self.npins if self.mirrorinputs else self.ppins

    @property
    def path(self):
        return ((-0.5, -0.5), (-0.5, 0.5), (0.5, 0))

    def draw(self, **kwargs):

        s = super (Uopamp, self).draw(**kwargs)

        if not self.nolabels:
            if self.mirrorinputs:
                s += self.annotate(self.node('lin+').s, '$-$', bold=True)
                s += self.annotate(self.node('lin-').s, '$+$', bold=True)
            else:
                s += self.annotate(self.node('lin+').s, '$+$', bold=True)
                s += self.annotate(self.node('lin-').s, '$-$', bold=True)
        return s    

    
class Ufdopamp(Chip):
    """This is for a fully differential opamp created with the U netlist
    type.  It has no wires.  See also Efdopamp for a fully differential
    opamp created with the E netlist type.

    """
    can_mirror = True

    auxiliary = {'lin+' : ('c', -0.375, 0.2),
                 'lin-' : ('c', -0.375, -0.2),
                 'lout+' : ('c', 0, -0.17),
                 'lout-' : ('c', 0, 0.17)}
                 
    auxiliary.update(Chip.auxiliary)
    
    ppins = {'out-' : ('r', 0.1, 0.2),
             'out+' : ('r', 0.1, -0.2),             
             'in+' : ('l', -0.5, 0.2),
             'in-' : ('l', -0.5, -0.2),
             'vdd' : ('t', -0.1, 0.3),
             'vss' : ('b', -0.1, -0.3),
             'vocm' : ('l', -0.5, 0)}

    npins = {'out+' : ('r', 0.1, +0.2),
             'out-' : ('r', 0.1, -0.2),             
             'in-' : ('l', -0.5, 0.2),
             'in+' : ('l', -0.5, -0.2),
             'vdd' : ('t', -0.1, 0.3),
             'vss' : ('b', -0.1, -0.3),
             'vocm' : ('l', -0.5, 0)}    

    pinlabels = {'vdd' : 'VDD', 'vss' : 'VSS'}

    @property
    def pins(self):
        return self.npins if self.mirrorinputs else self.ppins

    @property
    def path(self):
        return ((-0.5, -0.5), (-0.5, 0.5), (0.5, 0))

    def draw(self, **kwargs):

        s = super (Ufdopamp, self).draw(**kwargs)

        if not self.nolabels:
            if self.mirrorinputs:            
                s += self.annotate(self.node('lin+').s, '$-$', bold=True)
                s += self.annotate(self.node('lin-').s, '$+$', bold=True)
                s += self.annotate(self.node('lout+').s, '$-$', bold=True)
                s += self.annotate(self.node('lout-').s, '$+$', bold=True)
            else:
                s += self.annotate(self.node('lin+').s, '$+$', bold=True)
                s += self.annotate(self.node('lin-').s, '$-$', bold=True)
                s += self.annotate(self.node('lout+').s, '$+$', bold=True)
                s += self.annotate(self.node('lout-').s, '$-$', bold=True)
        return s

    
class Uinamp(Uopamp):
    """Instrumentation amplifier created with U netlist type."""

    can_mirror = True

    auxiliary = {'lin+' : ('c', -0.375, 0.3),                        
                 'lin-' : ('c', -0.375, -0.3)}
    auxiliary.update(Chip.auxiliary)
    
    ppins = {'out' : ('r', 0.5, 0.0),
             'in+' : ('l', -0.5, 0.3),
             'in-' : ('l', -0.5, -0.3),
             'vdd' : ('t', 0, 0.25),
             'vdd2' : ('t', -0.225, 0.365),
             'vss2' : ('b', -0.225, -0.365),
             'vss' : ('b', 0, -0.25),
             'ref' : ('b', 0.225, -0.135),
             'r+' : ('l', -0.5, 0.2),
             'r-' : ('l', -0.5, -0.2)}

    npins = {'out' : ('r', 0.5, 0.0),
             'in-' : ('l', -0.5, 0.3),
             'in+' : ('l', -0.5, -0.3),
             'vdd' : ('t', 0, 0.25),
             'vdd2' : ('t', -0.225, 0.365),
             'vss2' : ('b', -0.225, -0.365),
             'vss' : ('b', 0, -0.25),
             'ref' : ('b', 0.225, -0.135),
             'r-' : ('l', -0.5, 0.2),
             'r+' : ('l', -0.5, -0.2)}    


class Uisoamp(Ufdopamp):
    """Isolated amplifier created with U netlist type."""

    auxiliary = {'lin+' : ('c', -0.4, 0.2),
                 'lin-' : ('c', -0.4, -0.2),
                 'lout+' : ('c', 0.1, 0.12),
                 'lout-' : ('c', 0.1, -0.12)}    
    auxiliary.update(Chip.auxiliary)
    
    pins = {'out-' : ('r', 0.2, -0.15),
            'out+' : ('r', 0.2, 0.15),
            'out' : ('r', 0.5, 0),                         
            'in+' : ('l', -0.5, 0.2),
            'in-' : ('l', -0.5, -0.2),
            'in' : ('l', -0.5, 0.0),            
            'vdd1' : ('t', -0.3, 0.4),
            'vss1' : ('b', -0.3, -0.4),
            'vdd2' : ('t', 0.0, 0.25),
            'vss2' : ('b', 0.0, -0.25)}

    @property
    def path(self):
        return [((-0.5, -0.5), (-0.5, 0.5), (-0.2, 0.35), (-0.2, -0.35)),
                ((-0.1, -0.3), (-0.1, 0.3), (0.5, 0), (-0.1, -0.3))]
    
class Wire(Bipole):

    def __init__(self, sch, namespace, defname, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        opts = Opts(opts_string)
        
        implicit = False
        for key in self.implicit_keys:
            if key in opts:
                implicit = True
                break

        if implicit:
            # Rename second node since this is spatially different from
            # other nodes of the same name.  Add underscore at start so node
            # not drawn.

            # For example, the following nets share the ground node 0:
            # W 1 0; implicit
            # W 2 0; implicit            
            
            node_names = (node_names[0], '_' + name + '@' + node_names[1])
        super(Wire, self).__init__(sch, namespace, defname, name, cpt_type,
                                   cpt_id, string,
                                   opts_string, node_names, keyword, *args)
        self.implicit = implicit

    def draw_implicit(self, **kwargs):
        """Draw implicit wires, i.e., connections to ground, etc."""

        kind = None
        for key in self.implicit_keys:
            if key in self.opts:
                kind = key
                break;

        # I like the sground symbol for power supplies but rground symbol
        # is also common.
        if (kind is None) or (kind == 'implicit'):
            kind = 'sground'
        anchor = 'south west'
        if self.down:
            anchor = 'north west'

        n1, n2 = self.nodes
        s = self.draw_path((n1.s, n2.s))

        label_pos = n2
        if kind in ('input', 'output', 'bidir', 'pad'):

            # TODO: create node shapes
            
            try:
                scale = float(self.opts[kind])
            except:
                scale = 1.0

            h = 0.5 * scale * 1.2
            w = 0.75 * scale * 1.2
            a = 0.25 * scale * 1.2

            x = (w + a) / 2
            
            if kind == 'output':
                q = self.tf(n2.pos, ((0, h / 2), (w, h / 2),
                                     (w + a, 0), (w, -h / 2), (0, -h / 2)),
                            scale=1)
            elif kind == 'input':
                q = self.tf(n2.pos, ((0, 0), (a, h / 2), (a + w, h / 2),
                                     (a + w, -h / 2), (a, -h / 2)),
                            scale=1)
            elif kind == 'bidir':
                q = self.tf(n2.pos, ((0, 0), (a, h / 2), (w, h / 2),
                                     (a + w, 0), (w, -h / 2),
                                     (a, -h / 2)),
                            scale=1)                            
            elif kind == 'pad':
                q = self.tf(n2.pos, ((0, h / 2), (w + a, h / 2),
                                     (w + a, -h / 2), (0, -h / 2)),
                            scale=1)
            
            s += self.draw_path(q, closed=True)

            if 'l' in self.opts:
                lpos = self.tf(n2.pos, (x, 0), scale=1)
                s += r'  \draw[align=center] (%s) node {%s};''\n' % (
                    lpos, self.label(**kwargs))                
            return s

        rotate = 90
        if kind in ('antenna', 'rxantenna', 'txantenna'):
            rotate = 0

        args_str = 'scale=0.5, rotate=%d' % (self.angle + rotate)
        # mirror and invert keywords are ignored by circuitikz
        if self.mirror:
            args_str += ', yscale=-1'
        if self.invert:
            args_str += ', xscale=-1'            
        
        s += r'  \draw (%s) node[%s, %s] {};''\n' % (n2.s, kind, args_str)

        if 'l' in self.opts:
            lpos = self.tf(n2.pos, (0.125, 0), scale=1)
            s += r'  \draw[anchor=%s] (%s) node {%s};''\n' % (
                anchor, lpos, self.label(**kwargs))
        return s

    def setup(self):
        super(Wire, self).setup()

        if self.implicit:
            self.nodes[1].implicit = True        
    
    def draw(self, **kwargs):

        if not self.check():
            return ''

        if self.implicit:
            return self.draw_implicit(**kwargs)
            
        def arrow_map(name):

            try:
                return {'tee': '|', 'otri' : 'open triangle 60',
                        'tri' : 'triangle 60'}[name]
            except:
                return name

        n1, n2 = self.nodes

        # W 1 2; up, arrow=tri, l=V_{dd}
        # W 1 3; right, arrow=otri
        # W 1 4; down, arrow=tee, l=0V
        # W 1 5; left, startarrow=tri, endarrow=open triangle 90, bus=8

        startarrow = self.opts.pop('startarrow', '')
        endarrow = self.opts.pop('arrow', '')
        endarrow = self.opts.pop('endarrow', endarrow)

        bus = self.opts.pop('bus', False)
        style = ''
        if bus:
            # TODO if bus has numeric arg, indicate number of lines with slash.
            style = 'ultra thick'

        # TODO, add arrow shapes for earth symbol.

        if self.steps is None:
            s = r'  \draw[%s-%s, %s, %s] (%s) to (%s);''\n' % (
                arrow_map(startarrow), arrow_map(endarrow), style,
                self.args_str, n1.s, n2.s)
        else:

            steps = Steps(self.steps, n1.pos, n2.pos)

            path = '(%s)' % n1.pos
            for pos in steps:
                path += ' to (%s)' % pos

            s = r'  \draw[%s-%s, %s, %s] %s;''\n' % (
                arrow_map(startarrow), arrow_map(endarrow), style,
                self.args_str, path)

        if self.voltage_str != '':
            print('There is no voltage drop across an ideal wire!')

        if self.current_str != '' or self.label_str != '':
            # FIXME, we don't really want the wire drawn since this
            # can clobber the arrow.  We just want the current
            # annotation and/or the label.

            # To handle multiple labels, we need to draw separate wires.
            for label_str in self.label_str_list:
                s += r'  \draw[%s] (%s) [short, %s, %s] to (%s);''\n' % (
                    self.args_str, n1.s, self.current_str, label_str, n2.s)
            if self.label_str_list == []:
                s += r'  \draw[%s] (%s) [short, %s] to (%s);''\n' % (
                    self.args_str, n1.s, self.current_str, n2.s)
        return s


class FB(Bipole):
    """Ferrite bead"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes
        centre = (n1.pos + n2.pos) * 0.5
        w = 0.125
        h = 0.2
        
        q1 = self.tf(centre, ((-0.5 * w, -h), (-0.5 * w, h),
                              (0.5 * w, h), (0.5 * w, -h)), -30)
        q = self.tf(centre, ((-0.53 * w, 0), (0.53 * w, 0), (-w, -2 * h)))

        s = self.draw_path(q1, closed=True, style='thick')
        s += self.draw_path((n1.s, q[0]))
        s += self.draw_path((q[1], n2.s))
        s += self.draw_label(q[2], **kwargs)
        return s


class CPE(Bipole):
    """Constant phase element"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes
        centre = (n1.pos + n2.pos) * 0.5
        w = 0.125
        h = 0.2
        
        q1 = self.tf(centre, ((0, -h), (-w, 0), (0, h)), 0)
        q2 = self.tf(centre, ((w, -h), (0, 0), (w, h)), 0)
        q = self.tf(centre, ((-w, 0), (0, 0), (-w, -2 * h), (-w, 2 * h)))

        s = self.draw_path(q1, closed=False, style='thick')
        s += self.draw_path(q2, closed=False, style='thick')
        s += self.draw_path((n1.s, q[0]))
        s += self.draw_path((q[1], n2.s))
        # FOO
        s += self.draw_label(q[2], ('l', 'l_'), default=True, **kwargs)
        s += self.draw_label(q[3], ('l^', ), default=False, **kwargs)        
        return s    


class XX(Cpt):
    directive = True

    def draw(self, **kwargs):

        if self.string.startswith(';;'):
            return ' ' + self.string[2:] + '\n'
        return ''

    
classes = {}

def defcpt(name, base, docstring, cpt=None):
    
    if isinstance(base, str):
        base = classes[base]

    newclass = type(name, (base, ), {'__doc__': docstring})

    if cpt is not None:
        newclass.tikz_cpt = cpt
    classes[name] = newclass


def make(classname, parent, name, cpt_type, cpt_id,
         string, opts_string, node_names, *args):

    # Create instance of component object
    try:
        newclass = getattr(module, classname)
    except:
        newclass = classes[classname]

    cpt = newclass(parent, name, cpt_type, cpt_id, string, opts_string, 
                   node_names, *args)
    # Add named attributes for the args?   Lname1, etc.
        
    return cpt

# Dynamically create classes.
defcpt('ADC', Bipole, 'ADC', 'adc')
defcpt('AM', Bipole, 'Ammeter', 'ammeter')

defcpt('BAT', Bipole, 'Battery', 'battery')

defcpt('C', Bipole, 'Capacitor', 'C')

defcpt('D', Bipole, 'Diode', 'D')
defcpt('DAC', Bipole, 'DAC', 'dac')
defcpt('Dled', 'D', 'LED', 'leD')
defcpt('Dphoto', 'D', 'Photo diode', 'pD')
defcpt('Dschottky', 'D', 'Schottky diode', 'zD')
defcpt('Dtunnel', 'D', 'Tunnel diode', 'tD')
defcpt('Dzener', 'D', 'Zener diode', 'zD')

defcpt('E', VCS, 'VCVS', 'american controlled voltage source')
defcpt('F', VCS, 'CCCS', 'american controlled current source')
defcpt('G', CCS, 'VCCS', 'american controlled current source')
defcpt('H', CCS, 'CCVS', 'american controlled voltage source')

defcpt('FS', Bipole, 'Fuse', 'fuse')

defcpt('GY', Gyrator, 'Gyrator', 'gyrator')

defcpt('I', Bipole, 'Current source', 'I')
defcpt('sI', Bipole, 'Current source', 'I')
defcpt('Isin', 'I', 'Sinusoidal current source', 'sI')
defcpt('Idc', 'I', 'DC current source', 'I')
defcpt('Istep', 'I', 'Step current source', 'I')
defcpt('Iac', 'I', 'AC current source', 'sI')
defcpt('Inoise', 'I', 'Noise current source', 'sI')

defcpt('J', JFET, 'N JFET transistor', 'njfet')
defcpt('Jnjf', 'J', 'N JFET transistor', 'njfet')
defcpt('Jpjf', 'J', 'P JFET transistor', 'pjfet')

defcpt('L', Bipole, 'Inductor', 'L')

defcpt('M', MOSFET, 'N MOSJFET transistor', 'nmos')
defcpt('Mnmos', 'M', 'N channel MOSJFET transistor', 'nmos')
defcpt('Mpmos', 'M', 'P channel MOSJFET transistor', 'pmos')
defcpt('MISC', Bipole, 'Generic circuitikz component', '')

defcpt('NR', Bipole, 'Noiseless resistor', 'R')

defcpt('O', Bipole, 'Open circuit', 'open')
defcpt('P', Bipole, 'Port', 'open')

defcpt('Q', BJT, 'NPN transistor', 'npn')
defcpt('Qpnp', 'Q', 'PNP transistor', 'pnp')
defcpt('Qnpn', 'Q', 'NPN transistor', 'npn')

defcpt('R', Bipole, 'Resistor', 'R')
defcpt('RV', Potentiometer, 'Potentiometer', 'pR')

defcpt('Sbox', Box, 'Box shape')
defcpt('Scircle', Circle, 'Circle shape')
defcpt('Sellipse', Ellipse, 'Ellipse shape')
defcpt('Striangle', Triangle, 'Triangle shape')

defcpt('SW', Bipole, 'Switch', 'closing switch')
defcpt('SWno', 'SW', 'Normally open switch', 'closing switch')
defcpt('SWnc', 'SW', 'Normally closed switch', 'opening switch')
defcpt('SWpush', 'SW', 'Pushbutton switch', 'push button')
defcpt('SWspdt', SPDT, 'SPDT switch', 'spdt')

defcpt('TF', Transformer, 'Transformer', 'ideal transformer')
defcpt('TFcore', Transformer, 'Transformer with core', 'transformer core')
defcpt('TFtapcore', TFtap, 'Tapped transformer with core', 'transformer core')
defcpt('TP', TwoPort, 'Two port', '')

defcpt('Ubox', Box2, 'Box')
defcpt('Ucircle', Circle2, 'Circle')
defcpt('Ubox4', Box4, 'Box')
defcpt('Ubox12', Box12, 'Box')
defcpt('Ucircle4', Circle4, 'Circle')

defcpt('V', Bipole, 'Voltage source', 'V')
defcpt('sV', Bipole, 'Voltage source', 'V')
defcpt('Vsin', 'V', 'Sinusoidal voltage source', 'sV')
defcpt('Vdc', 'V', 'DC voltage source', 'V')
defcpt('Vstep', 'V', 'Step voltage source', 'V')
defcpt('Vac', 'V', 'AC voltage source', 'sV')
defcpt('Vnoise', 'V', 'Noise voltage source', 'sV')

defcpt('VCVS', VCS, 'VCVS', 'american controlled voltage source')
defcpt('CCCS', VCS, 'CCCS', 'american controlled current source')
defcpt('VCCS', CCS, 'VCCS', 'american controlled current source')
defcpt('CCVS', CCS, 'CCVS', 'american controlled voltage source')

defcpt('VM', Bipole, 'Voltmeter', 'voltmeter')

defcpt('W', Wire, 'Wire', 'short')

defcpt('XT', Bipole, 'Crystal', 'piezoelectric')

defcpt('Y', Bipole, 'Admittance', 'european resistor')
defcpt('Z', Bipole, 'Impedance', 'european resistor')

defcpt('m', Bipole, 'Mass', 'mass')
defcpt('k', Bipole, 'Spring', 'spring')
defcpt('r', Bipole, 'Damper', 'damper')

# Perhaps AM for ammeter, VM for voltmeter, VR for variable resistor?
# Currently, a variable resistor is supported with the variable
# option.
