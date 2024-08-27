from .fixedcpt import FixedCpt


class Shape(FixedCpt):
    """General purpose shape"""

    default_aspect = 1.0
    can_mirror = True
    can_invert = True
    pinlabels = {}

    auxiliary = {'mid': ('c', 0.0, 0.0),
                 'bl': ('l', -0.5, -0.5),
                 'br': ('r', 0.5, -0.5),
                 'top': ('t', 0, 0.5),
                 'tl': ('l', -0.5, 0.5),
                 'tr': ('r', 0.5, 0.5)}

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
        # pinlabels=default label the default pinlabels
        # pinlabels=none label no pinlabels

        def pinlabel(node_name):

            fields = node_name.split('.')
            pinname = fields[-1]

            try:
                return self.pinlabels[pinname]
            except:
                pass
            return ''

        prefix = self.name + '.'

        pinlabels = self.opts.get('pinlabels', 'default')

        if pinlabels == 'none':
            return {}
        elif pinlabels in ('', 'auto', 'connected'):
            return {name: pinlabel(name) for name in self.ref_node_names if pinlabel(name) != ''}
        elif pinlabels == 'all':
            return {name: pinlabel(name) for name in self.pin_node_names if pinlabel(name) != ''}
        elif pinlabels == 'default':
            return {self.name + '.' + name: pinlabel(self.name + '.' + name) for name in self.default_pins if pinlabel(name) != ''}
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

        for node_name, pinlabel in pinlabels.items():
            # Add pin to nodes so that it will get allocated a coord.
            node = self.sch._node_add(node_name, self, auxiliary=True)
            node.pin = not self.fakepin(node.basename)
            node.pinpos = self.pinpos(node.basename)

            # TODO, perhaps use pinlabel to indicate clock?
            node.clock = pinlabel != '' and pinlabel[0] == '>'
            if node.clock:
                # Remove clock designator
                pinlabel = pinlabel[1:]

            node.pinlabel = pinlabel

    def process_pinnodes(self):

        self.process_nodes(self.parse_pinnodes(), draw_pin=True)

    def process_pinnames(self):

        self.process_nodes(self.parse_pinnames(), add_pinname=True)

    def process_implicit_nodes(self):

        for pin_name in self.pin_node_names:
            if pin_name not in self.sch.nodes:
                continue
            node = self.sch.nodes[pin_name]
            implicit = self.implicit_key(node.opts)
            if implicit:
                node.implicit = True
                node.implicit_symbol = implicit
                node.pinpos = self.pinpos(node.basename)

    def setup(self):

        super(Shape, self).setup()

        self.process_pinnodes()
        self.process_pinlabels()
        self.process_pinnames()

        # Ensure all the shape nodes are marked as pins.
        for node in self.nodes:
            node.pin = not self.fakepin(node.basename)

    def draw(self, scale=1, **kwargs):

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

        cpt_args_str = self.cpt_args_str(**kwargs)
        if not self.nodraw:
            cpt_args_str += ', draw'

        # shape border rotate rotates the box but not the text
        s = r'  \draw (%s) node[%s, thick, inner sep=0pt, minimum width=%.2fcm, minimum height=%.2fcm, text width=%.2fcm, align=center, shape border rotate=%s, %s, scale=%s] (%s) {%s};''\n' % (
            self.centre, self.shape, self.width, self.height,
            text_width, self.angle, cpt_args_str, scale, self.s, label)
        return s
