from warnings import warn
from numpy import zeros, array, pi, cos, sin, dot
from ...label import Label
from ...labels import Labels
from ...opts import Opts
from ...parser import split
from ...schemmisc import Pos, Steps
from ...label import Label
from ...labels import Labels
from ...latex import latex_format_label
from ...config import implicit_default
from ..utils import check_boolean


def anchor_map(anchor):

    mapping = {'n': 'north', 'e': 'east',
               's': 'south', 'w': 'west',
               'ne': 'north east', 'nw': 'north west',
               'se': 'south east', 'sw': 'south west'}

    try:
        anchor = mapping[anchor]
    except KeyError:
        pass
    return anchor


def anchor_choose(pinpos, outside=False):

    if outside:
        anchors = {None: 'south east',
                   'c': 'south east',
                   'l': 'west', 'r': 'east',
                   't': 'north', 'b': 'south'}
    else:
        anchors = {None: 'south east',
                   'c': 'south east',
                   'l': 'south east', 'r': 'north west',
                   't': 'south west', 'b': 'north west'}

    return anchors[pinpos]


def angle_choose(pinpos):

    angle = {'l': 180, 't': 90, 'b': -90, 'r': 0, None: 0}[pinpos]
    return angle


class Cpt(object):

    voltage_keys = ('v', 'v_', 'v^', 'v_>', 'v_<', 'v^>', 'v^<',
                    'v<', 'v>')
    current_keys = ('i', 'i_', 'i^', 'i_>',  'i_<', 'i^>', 'i^<',
                    'i>_', 'i<_', 'i>^', 'i<^', 'i>', 'i<', 'ir')
    flow_keys = ('f', 'f_', 'f^', 'f_>',  'f_<', 'f^>', 'f^<',
                 'f>_', 'f<_', 'f>^', 'f<^', 'f>', 'f<')
    label_keys = ('l', 'l_', 'l^')
    label2_keys = ('l2', 'l2_', 'l2^')
    annotation_keys = ('a', 'a_', 'a^')
    inner_label_keys = ('t', )
    connection_keys = ('input', 'output', 'bidir', 'pad')
    ground_keys = ('ground', 'sground', 'rground',
                   'cground', 'nground', 'pground', '0V')
    supply_positive_keys = ('vcc', 'vdd')
    supply_negative_keys = ('vee', 'vss')
    supply_keys = supply_positive_keys + supply_negative_keys
    implicit_keys = ('implicit', ) + ground_keys + supply_keys
    # The following keys do not get passed through to circuitikz.
    # Note, thickness is currently ignored.
    misc_keys = ('left', 'right', 'up', 'down', 'rotate', 'size',
                 'mirror', 'invert', 'scale', 'invisible', 'variable', 'fixed',
                 'aspect', 'pins', 'image', 'offset', 'pinlabels',
                 'pinnames', 'pinnodes', 'pindefs', 'outside',
                 'pinmap', 'kind', 'wire', 'ignore', 'style', 'nosim',
                 'nowires', 'nolabels', 'steps', 'free', 'fliplr', 'flipud',
                 'nodots', 'draw_nodes', 'label_nodes', 'nodraw',
                 'mirrorinputs', 'autoground', 'xoffset', 'yoffset',
                 'anchor', 'def', 'nodes', 'shape')
    label_opt_keys = ('label_values', 'label_ids', 'annotate_values',
                      'label_style', 'label_flip')
    color_keys = ('blue', 'cyan', 'magenta', 'yellow',
                  'black', 'gray', 'white', 'darkgray',
                  'lightgray', 'brown', 'lime', 'olive',
                  'orange', 'pink', 'purple', 'teal', 'violet')
    linestyle_keys = ('ultra thin', 'very thin', 'thin', 'semithick',
                      'thick', 'very thick', 'ultra thick',
                      'solid', 'dotted', 'densely dotted',
                      'loosely dotted', 'dashed', 'densely dashed',
                      'loosely dashed')

    all_label_keys = voltage_keys + current_keys + flow_keys + \
        label_keys + label2_keys + inner_label_keys + \
        annotation_keys

    special_keys = all_label_keys + misc_keys + \
        implicit_keys + label_opt_keys + connection_keys

    node_special_keys = ('label', 'l', 'anchor')

    can_rotate = True
    can_scale = False
    can_mirror = False
    can_invert = False
    do_transpose = False
    can_stretch = True
    shape_scale = 1.0
    default_width = 1.0
    default_aspect = 1.0
    # node_pinnames maps node numbers to pinnames
    node_pinnames = ()
    default_pins = ()
    pins = {}
    # Auxiliary nodes are used for finding the centre of the shape or
    # to define a bounding box.
    auxiliary = {}
    required_auxiliary = ('mid', )
    directive = False
    place = True
    kinds = {}
    styles = {}
    aliases = {}

    @property
    def s(self):
        """Sanitised name"""
        return self.name.replace('.', '@')

    def __init__(self, sch, namespace, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        self.sch = sch
        self.type = cpt_type
        self.id = cpt_id
        self.name = name
        self.namespace = namespace

        self.string = string
        self.net = string.split(';')[0]
        self.opts_string = opts_string
        self.opts = Opts(self.opts_string)
        self._setup = False

        self.args = args
        self.classname = self.__class__.__name__

        # Drawing hints
        self.opts = Opts(opts_string)

        self._process_opts()

        self.labels = Labels()
        for key, val in self.opts.items():
            if key in self.all_label_keys:
                self.labels.add(key, val)

        # List of node names specified as component arguments.  For example,
        # given R 1 2, node_names = ['1', '2'].
        # The ordering of this list is important.
        self.node_names = list(node_names)

        prefix = self.name + '.'

        auxiliary_node_names = []
        for pin in self.auxiliary:
            auxiliary_node_names.append(prefix + pin)

        self.auxiliary_node_names = auxiliary_node_names


    def _process_opts(self):

        if self.sch is None:
            return

        defines = self.sch.defines

        # Look for defs
        if 'def' in self.opts:
            for arg in self.opts['def']:
                parts = arg.split('=', 1)
                arg = parts[1]
                if arg.startswith('{') and arg.endswith('}'):
                    arg = arg[1:-1]
                defines[parts[0]] = arg

        # Replace defs
        for opt in self.opts.copy():
            if opt in defines:
                self.opts.pop(opt)

                defs = split(defines[opt], ',')

                for def1 in defs:
                    def1 = def1.strip()
                    parts = def1.split('=', 1)
                    if len(parts) == 2:
                        self.opts[parts[0]] = parts[1]
                    else:
                        self.opts[parts[0]] = ''

    def check_nodes(self):

        # There are 5 cases:
        # 1. node             R1 1 2
        # 2. pin ref          R1 1 U1.in
        # 3. include node     R1 1 s.2
        # 4. include pin ref  R1 1 s.U1.in
        # 5. relative ref     R1 1 ._2   or  R1 1 R1._2
        # Note R1 1 ._2 gets converted to R1 1 R1._2 before this is called.

        for node_name in self.node_names:

            fields = node_name.split('.')
            # Case 1
            if len(fields) < 2:
                continue

            cpt_name = '.'.join(fields[0:-1])
            basename = fields[-1]

            # Case 3
            if cpt_name in self.sch.subnetlists:
                continue

            if cpt_name not in self.sch.elements:
                raise ValueError('Undefined component `%s` in `%s`' %
                                 (cpt_name, self))
            cpt = self.sch.elements[cpt_name]

            # Case 5
            if not isinstance(cpt, Shape):
                continue

            if node_name not in cpt.all_node_names:
                known_pins = ', '.join(cpt.all_node_names)
                raise ValueError(
                    'Unknown node `%s` for `%s` in `%s`, known nodes: %s' %
                    (basename, cpt_name, self, known_pins))

        self.relative_node_names = []
        for name in self.node_names:
            fields = name.split('.')
            if len(fields) < 2:
                continue
            if fields[-2] == self.name:
                self.relative_node_names.append(name)

    def parse_nodes(self):

        self.allpins = self.pins.copy()
        self.allpins.update(self.auxiliary)

        prefix = self.name + '.'
        pin_node_names = []
        for pin in self.pins.keys():
            pin_node_names.append(prefix + pin)

        self.pindefs = self.parse_pindefs()
        for node_name, pindef in self.pindefs.items():
            if True:
                # Remove old name.
                self.allpins[pindef] = self.allpins.pop(node_name)
                pin_node_names.remove(prefix + node_name)
            else:
                # Keep old name as an alias.  This will cause problems
                # when printing pinnodes.
                self.allpins[pindef] = self.allpins[node_name]
            pin_node_names.append(prefix + pindef)

        # These are all the pin names belonging to the cpt.
        self.pin_node_names = pin_node_names

        # These are all the pin nodes required to be shown for the cpt.
        # This is set by the process_pins method.
        self.drawn_pins = []

        attribute_node_names = []
        if isinstance(self, Shape):
            for opt in self.node_opts():
                if opt[0] not in self.node_names:
                    attribute_node_names.append(prefix + opt[0])

        self.attribute_node_names = attribute_node_names

        alias_node_names = []
        for alias in self.aliases:
            alias_node_names.append(prefix + alias)

        self.all_node_names = self.required_node_names + \
            self.auxiliary_node_names + pin_node_names + \
            attribute_node_names + alias_node_names

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
        return float(val) * self.shape_scale

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

    @property
    def xoffset(self):
        return float(self.opts.get('xoffset', 0))

    @property
    def yoffset(self):
        return float(self.opts.get('yoffset', 0))

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
        return self.boolattr('invert') or self.fliplr

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
        """Normalized width"""
        return 1.0

    @property
    def h(self):
        """Normalized height"""
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

    def anchor_opt(self, thing, default=None):

        val = thing.opts.get('anchor', default)
        return anchor_map(val)

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
            return array(Rdict[angle])

        t = angle / 180.0 * pi
        return array(((cos(t), sin(t)),
                         (-sin(t), cos(t))))

    @property
    def required_node_names(self):
        """Subset of node_names.  This filters out nodes that are not
        drawn.  For example, the ground node of an Eopamp is not drawn."""

        # The node_pinnames tuple specifies the required nodes.
        node_names = []
        for pinname, node_name in zip(self.node_pinnames, self.node_names):
            if pinname != '':
                node_names.append(node_name)
        return node_names

    def node(self, pinname):
        """Return node by pinname"""

        if pinname in self.aliases:
            pinname = self.aliases[pinname]

        if pinname in self.node_pinnames:
            index = self.node_pinnames.index(pinname)
            node_name = self.node_names[index]
        else:
            node_name = self.name + '.' + pinname

        try:
            return self.sch.nodes[node_name]
        except KeyError:
            known_pins = ', '.join(self.pins.keys())
            raise ValueError('Unknown pinname %s for `%s`, known pins: %s' %
                             (pinname, self, known_pins))

    def pinpos(self, pinname):
        """Return pinpos by pinname"""

        if pinname in self.aliases:
            pinname = self.aliases[pinname]

        if pinname not in self.allpins:
            raise ValueError('Unknown pin %s for %s, known pins: %s'
                             % (pinname, self.name, ', '.join(list(self.allpins))))
        return self.allpins[pinname][0][0]

    def fakepin(self, pinname):
        """Return True if pinname is a fake pin"""

        if pinname in self.aliases:
            pinname = self.aliases[pinname]

        # Check if known pinname; it might be a netlist node name
        if pinname not in self.allpins:
            return True

        # Fake pins have a pinpos like lx, rx, tx, or bx
        return len(self.allpins[pinname][0]) == 2

    def node_opts(self):

        opts = []
        for key, val in self.opts.items():
            if not key.startswith('.'):
                if '.' in key:
                    warn('The node key %s should start with a .' % key)
                continue

            parts = key[1:].split('.')
            if len(parts) != 2:
                raise ValueError('Badly formatted node attribute: ' + key)
            node = parts[0]

            if node in self.aliases:
                node = self.aliases[node]

            if node not in self.pins:
                raise ValueError('Unknown pin %s in %s.  Known pins: %s' %
                                 (node, key, ', '.join(self.pins)))
            parts.append(val)
            opts.append(parts)
        return opts

    @property
    def nodes(self):
        """Nodes used to draw the element."""

        if not self._setup:
            raise RuntimeError('Using nodes before setup.')

        if hasattr(self, '_nodes'):
            return self._nodes

        # Perhaps determine coords here as well and cache them?

        # The ordering of the first nodes is important.
        # These must match with the electrical nodes.
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
                if pinname in self.aliases:
                    pinname = self.aliases[pinname]
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

        tcoords = zeros((len(coords), 2))
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

    def draw_args(self, opts, **kwargs):

        allow_keys = self.linestyle_keys

        return opts.as_list(allow=allow_keys, **kwargs)

    def cpt_args(self, opts, **kwargs):

        # Allow draw, fill, color

        ignore_keys = self.special_keys

        return opts.as_list(ignore_keys, **kwargs)

    def draw_cptnode(self, pos, cpt='', args='', dargs='', label=''):
        """Create a string to draw a tikz node containing a circuitikz
        component `cpt` at position `pos`. `args` is a list or string
        of the node options; `dargs` is a list or string of the draw
        options.  `label` is an optional label.

        The general form of the generated string is:

        \draw[dargs] (pos) node[cpt, args] {label};

        """

        if isinstance(args, list):
            args = ', '.join([arg for arg in args if arg != ''])
        if isinstance(dargs, list):
            dargs = ', '.join([arg for arg in dargs if arg != ''])
        label_str = latex_format_label(label)

        if args == '':
            args = cpt
        elif cpt != '':
            args = cpt + ', ' + args

        if args == '':
            node_str = 'node'
        else:
            node_str = 'node[' + args + ']'

        if dargs == '':
            s = r'  \draw (%s) %s {%s};''\n' % (
                pos, node_str, label_str)
        else:
            s = r'  \draw[%s] (%s) %s {%s};''\n' % (
                dargs, pos, node_str, label_str)
        return s

    def draw_cpt(self, pos1, pos2, cpt='', args='', dargs=''):
        """Create a string to draw a circuitikz component `cpt` between
        positions `pos1` and `pos2`. `args` is a list or string of the
        component options; `dargs` is a list or string of the draw options.

        The general form of the generated string is:

        \draw[dargs] (pos1) to [cpt, args] (pos2);

        """

        if isinstance(args, list):
            args = ', '.join([arg for arg in args if arg != ''])
        if isinstance(dargs, list):
            dargs = ', '.join([arg for arg in dargs if arg != ''])

        if args == '':
            args = cpt
        elif cpt != '':
            args = cpt + ', ' + args

        if dargs != '':
            dargs = '[' + dargs + ']'
        if args != '':
            args = ' [' + args + ']'

        s = r'  \draw%s (%s) to%s (%s);''\n' % (
            dargs, pos1, args, pos2)
        return s

    def draw_path(self, points, style='', join='--', closed=False, dargs=None):

        path = (' %s ' % join).join(['(%s)' % point for point in points])
        if closed:
            path += ' %s cycle' % join

        if dargs is None:
            dargs = self.draw_args(self.opts)
        if style != '':
            dargs.append(style)
        dargs = ', '.join([arg for arg in dargs if arg != ''])
        if dargs != '':
            dargs = '[' + dargs + ']'

        s = r'  \draw%s %s;''\n' % (dargs, path)
        return s

    def draw_connection(self, n, kind):
        """Draw connection and label."""

        args = self.draw_args(n.opts)

        try:
            scale = float(self.opts[kind])
        except:
            scale = 1.0

        pinpos = n.pinpos
        angle = angle_choose(pinpos)

        h = 0.5 * scale * 1.2
        w = 0.75 * scale * 1.2
        a = 0.25 * scale * 1.2
        x = (w + a) / 2

        if kind == 'output':
            q = self.tf(n.pos, ((0, h / 2), (w, h / 2),
                                (w + a, 0), (w, -h / 2), (0, -h / 2)),
                        scale=1, angle_offset=angle)
        elif kind == 'input':
            q = self.tf(n.pos, ((0, 0), (a, h / 2), (a + w, h / 2),
                                (a + w, -h / 2), (a, -h / 2)),
                        scale=1, angle_offset=angle)
        elif kind == 'bidir':
            q = self.tf(n.pos, ((0, 0), (a, h / 2), (w, h / 2),
                                (a + w, 0), (w, -h / 2),
                                (a, -h / 2)),
                        scale=1, angle_offset=angle)
        elif kind == 'pad':
            q = self.tf(n.pos, ((0, h / 2), (w + a, h / 2),
                                (w + a, -h / 2), (0, -h / 2)),
                        scale=1, angle_offset=angle)

        s = self.draw_path(q, closed=True, dargs=args)

        label = n.opts.get('l', n.opts.get('label', ''))
        if label != '':

            opts = n.opts.copy()
            opts.pop('fill', None)
            args = self.draw_args(opts)

            lpos = self.tf(n.pos, (x, 0), scale=1, angle_offset=angle)
            dargs = ['align=center']
            angle += self.angle
            if angle > 90:
                angle -= 180
            elif angle < -90:
                angle += 180

            if angle != 0:
                args = ['rotate=%s' % angle]
            s += self.draw_cptnode(lpos, args=args, dargs=dargs, label=label)

        return s

    def draw_implicit_rotate(self, n, kind, draw_nodes):
        """Draw implicit connection and label with rotation."""

        label = ''
        if kind == '0V':
            label = r'0\,\mathrm{V}'
            kind = implicit_default
        elif kind == 'implicit':
            kind = implicit_default

        args = self.draw_args(n.opts)
        pinpos = n.pinpos

        if self.type != 'A':
            angle = angle_choose(pinpos)
        else:
            angle = 0

        if kind in ('vcc', 'vdd'):
            # vcc/vee and vdd/vss are drawn in opposite directions
            # so use vss to be consistent.
            kind = 'vss'

        # The default drawing direction for a ground symbol
        # is down (rotate = 0).  Need to rotate by 90 degrees
        # so convert Lcapy's direction to Circuitikz's direction.
        rotate = angle + self.angle + 90
        if rotate != 0:
            args.append('rotate=%d' % rotate)

        s = self.draw_cptnode(n.s, kind, args)

        # vss and vdd labels are drawn in the correct place except
        # with rotation.  However, there is no provision for labels
        # for ground, sground, etc.  So we position them ourself.

        label = n.opts.get('l', n.opts.get('label', label))
        if label != '':
            lpos = self.tf(n.pos, (0.5, 0), scale=1, angle_offset=0)

            if abs(self.angle) <= 90:
                anchor = 'west'
                xoffset = 0.2
            else:
                xoffset = -0.2
                anchor = 'east'
            lpos.x += xoffset

            # This does not rotate the text (unlike connections).
            dargs = ['anchor=' + anchor]
            s += self.draw_cptnode(lpos, dargs=dargs, label=label)

        return s

    def draw_implicit_norotate(self, n, kind, draw_nodes):
        """Draw implicit connection and label with no rotation."""

        label = ''
        if kind == '0V':
            label = r'0\,\mathrm{V}'
            kind = implicit_default
        elif kind == 'implicit':
            kind = implicit_default

        args = self.draw_args(n.opts)
        label = n.opts.get('l', n.opts.get('label', label))

        # vss and vdd labels are drawn in the correct place except
        # with rotation.  However, there is no provision for labels
        # for ground, sground, etc.  So we position them ourself.
        if kind in self.supply_keys:
            s = self.draw_cptnode(n.s, kind, args=args, label=label)
        else:
            s = self.draw_cptnode(n.s, kind, args=args)
            if label != '':

                lpos = self.tf(n.pos, (1, 0), scale=1,
                               angle_offset=-90)

                # This does not rotate the text (unlike connections).
                # TODO, use anchor instead of align
                dargs = ['align=center']
                s += self.draw_cptnode(lpos, dargs=dargs, label=label)

        return s

    def draw_implicit(self, n, kind, draw_nodes):
        """Draw implicit connection and label."""

        if kind in self.supply_keys:
            s = self.draw_implicit_norotate(n, kind, draw_nodes)
        else:
            s = self.draw_implicit_rotate(n, kind, draw_nodes)

        return s

    def draw_node(self, n, draw_nodes, dargs):
        """Draw a node symbol.  This also draws the node label
        for implicit and connection nodes."""

        # Give preference to annotation elts.
        if self.type != 'A':
            for elt in n.elt_list:
                if elt.type == 'A':
                    return ''

        if n.drawn:
            return ''
        n.drawn = True

        # Don't draw nodes for open-circuits.  Use port if want nodes drawn.
        if self.type == 'O':
            return ''

        args = self.draw_args(n.opts)

        s = ''

        kind = n.implicit_symbol
        if kind in self.connection_keys:
            s += self.draw_connection(n, kind)
        elif kind is not None:
            s += self.draw_implicit(n, kind, draw_nodes)

        if not n.visible(draw_nodes) or n.pin or not draw_nodes:
            return s

        symbol = n.opts.get('symbol', 'ocirc' if n.is_port or
                            n.is_dangling else 'circ')

        s += self.draw_cptnode(n.s, symbol, args, dargs)
        return s

    def draw_nodes(self, **kwargs):

        draw_nodes = self.draw_nodes_opt
        if draw_nodes is None:
            draw_nodes = kwargs.get('draw_nodes', True)
        dargs = self.draw_args(self.opts, **kwargs)

        s = ''
        for n in self.drawn_nodes:
            s += self.draw_node(n, draw_nodes=draw_nodes, dargs=dargs)

        return s

    def draw_pins(self):

        s = ''
        for n in self.drawn_pins:
            s += self.draw_cptnode(n.s, 'ocirc')
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

        anchor = anchor_choose(pinpos, outside != '')

        s = self.draw_cptnode(node.s, dargs='anchor=' + anchor,
                              label=node.pinlabel.replace('_', r'\_'))
        return s

    def draw_pinname(self, node):

        if node.pinname == '':
            return ''

        pinpos = node.pinpos
        pinpos = self.pinpos_rotate(pinpos, self.angle)

        # Move pinnames to outside
        mapping = {None: None,
                   'c': 'c',
                   'l': 'r', 'r': 'l',
                   't': 'b', 'b': 't'}
        pinpos = mapping[pinpos]

        anchor = anchor_choose(pinpos, True)

        s = self.draw_cptnode(node.s, dargs='anchor=' + anchor,
                              label=node.pinname.replace('_', r'\_'))
        return s

    def draw_node_label(self, node, label_nodes, anchor, dargs=None):

        if node.label_drawn:
            return ''
        node.label_drawn = True

        if not node.show_label(label_nodes):
            return ''

        anchor = self.anchor_opt(self, anchor)

        dargs = [] if dargs is None else dargs
        dargs.append('anchor=' + anchor)

        s = self.draw_cptnode(node.s, dargs=dargs, label=node.label)
        return s

    def draw_node_labels(self, **kwargs):

        label_nodes = self.label_nodes_opt
        if label_nodes is None:
            label_nodes = kwargs.get('label_nodes', 'primary')

        anchor = self.anchor_opt(self, kwargs.get('anchor', 'south east'))

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
                dargs = self.draw_args(self.opts, **kwargs)
                s += self.draw_node_label(node, label_nodes, anchor, dargs)

        return s

    def draw(self, **kwargs):
        raise NotImplementedError('draw method not implemented for %s' % self)

    def find_ref_node_names(self):
        """Determine which nodes are referenced for this component."""

        ref_node_names = []

        for node_name, node in self.sch.nodes.items():
            if node_name in self.relative_node_names:
                pass
            elif node.belongs(self.name):
                if node.basename not in self.allpins:
                    continue
                node.pin = not self.fakepin(node.basename)
                if node.basename not in self.auxiliary:
                    ref_node_names.append(node.name)
            elif self.namespace != '' and node_name.startswith(self.namespace):
                # Need to be lenient here since can have any old name.
                pass

        return ref_node_names

    def autoground(self, autoground):

        if (autoground not in (None, 'none') and
                autoground not in self.ground_keys):
            raise ValueError('Invalid autoground % s.  Choices are % s' %
                             (autoground, ', '.join(self.ground_keys)))

        for m, node_name in enumerate(self.required_node_names):
            if node_name != '0':
                continue

            new_node = self.sch.nodes[node_name].split(self)
            new_node.implicit = True
            new_node.implicit_symbol = autoground

            self.node_names[m] = new_node.name
            self.sch.nodes[new_node.name] = new_node

            index = self.all_node_names.index(node_name)
            self.all_node_names[index] = new_node.name

    def process_nodes(self, nodes, draw_pin=False, add_pinname=False):

        for n in nodes:
            # Add pin to nodes so that it will get allocated a coord.
            node = self.sch._node_add(n, self, auxiliary=True)

            node.pin = not self.fakepin(node.basename)
            node.pinpos = self.pinpos(node.basename)
            if draw_pin:
                self.drawn_pins.append(node)
            if add_pinname:
                node.pinname = node.basename

    def process_attribute_nodes(self):
        """Process nodes specified as attributes.  For example,
        R1 1 2; .n.vss"""

        self.process_nodes(self.attribute_node_names)

        opts = self.node_opts()
        for opt in opts:
            n = opt[0]
            node = self.node(n)
            node.opts.add(opt[1] + '=' + opt[2])

    def setup(self):
        self.ref_node_names = self.find_ref_node_names()
        self.process_attribute_nodes()

        self._setup = True

    def implicit_key(self, opts):

        prevkey = None
        for key in self.implicit_keys + self.connection_keys:
            if key in opts:
                if prevkey is not None:
                    raise ValueError(
                        'Multiple implicit node options %s and %s for %s: ' %
                        (prevkey, key, self))
                prevkey = key

        return prevkey

    def process_implicit_nodes(self):
        """Parse implicit nodes."""

        def split_nodes1(m, kind):
            node_name = self.node_pinnames[m]

            node = self.node(node_name)
            if node.pin:
                raise RuntimeError('Cannot split pin ' + node.name)

            if self.type != 'A':
                new_node = node.split(self)
            else:
                new_node = node

            new_node.implicit_symbol = kind
            new_node.implicit = True
            new_node.pinpos = self.pinpos(node_name)
            self.node_names[m] = new_node.name
            self.nodes[m] = new_node
            self.sch.nodes[new_node.name] = new_node

            index = self.all_node_names.index(node.name)
            self.all_node_names[index] = new_node.name
            return new_node

        # Old syntax, look for implicit, ground, etc.
        implicit = self.implicit_key(self.opts)
        if implicit:
            if implicit in self.supply_positive_keys:
                m = 0
            else:
                m = len(self.node_names) - 1
            new_node = split_nodes1(m, implicit)

            label = self.opts.pop('label', None)
            if label is not None:
                new_node.opts.add('l=' + label)
            label = self.opts.pop('l', None)
            if label is not None:
                new_node.opts.add('l=' + label)
            # Perhaps copy all the component opts to the node?
            for opt in ('anchor', 'fill', 'color'):
                val = self.opts.get(opt, None)
                if val is not None:
                    new_node.opts.add(opt + '=' + val)

        # New syntax, look for .p.implicit, .p.ground, etc.
        for m, node_name in enumerate(self.required_node_names):
            node = self.sch.nodes[node_name]
            implicit = self.implicit_key(node.opts)
            if implicit:
                split_nodes1(m, implicit)

    def parse_pindefs(self):
        return {}

    def draw_args_str(self, **kwargs):

        return ','.join(self.draw_args(self.opts, **kwargs))

    def cpt_args_str(self, **kwargs):

        return ','.join(self.cpt_args(self.opts, **kwargs))

    def label(self, keys=None, default=True, **kwargs):

        if keys is None:
            keys = self.label_keys

        label_values = check_boolean(kwargs.get('label_values', True))
        label_ids = check_boolean(kwargs.get('label_ids', True))

        label_str = ''
        if label_ids is True:
            label_str = self.id_label
        if label_values and self.value_label != '':
            label_str = latex_format_label(self.value_label)

        has_label = False
        for key, val in self.opts.items():
            if key in self.label_keys:
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

    def label_tweak(self, label, xscale, yscale, angle):

        # Circuitikz scales the label text so we undo this.
        if xscale != 1 or yscale != 1:
            label = r'\scalebox{%s}[%s]{%s}' % (1 / xscale, 1 / yscale, label)
        # Circuitikz rotates the label text so we undo this.
        if angle != 0:
            label = r'\rotatebox{%s}{%s}' % (-angle, label)
        return label

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
            raise ValueError(
                'Mutually exclusive drawing directions for %s' % self.name)

        return not self.invisible

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

        return centre + dot((x * self.w, y * self.h), self.R(angle_offset)) * scale

    def annotate(self, pos, label, dargs=None, bold=False):

        if bold:
            if label.startswith('$') and label.endswith('$'):
                # boldsymbol but requires amssym or amsmath.
                label = r'$\Large \boldsymbol{%s}$' % label[1:-1]
            else:
                label = r'\textbf{%s}' % label

        return self.draw_cptnode(pos, dargs=dargs, label=label)

    def draw_label(self, pos, keys=None, default=True, **kwargs):
        """Draw label for component that does not have a circuitikz label."""

        if keys is None:
            keys = self.label_keys

        args = self.draw_args(self.opts, **kwargs)

        return self.annotate(pos, self.label(keys, default=default, **kwargs),
                             args)


from .shape import Shape  # nopep8
