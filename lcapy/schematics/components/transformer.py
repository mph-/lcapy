from warnings import warn
from .fixedcpt import FixedCpt


class TF(FixedCpt):

    def draw_core(self, centre, l):

        s = ''

        core = 'two'
        if 'core' in self.opts:
            core = self.opts['core']

        if core in ('two', True):
            # Draw core with two lines
            q = self.tf(centre, ((-0.05, -l), (-0.05, l),
                                 (0.05, -l), (0.05, l)))
            s += self.draw_path(q[0:2], style='thick')
            s += self.draw_path(q[2:4], style='thick')
        elif core == 'one':
            # Draw core with one line
            q = self.tf(centre, ((0, -l), (0, l)))
            s += self.draw_path(q[0:2], style='thick')
        elif core not in ('', False):
            raise ValueError('Invalid core attribute value ' + core)
        return s


class TF1(TF):
    """Transformer"""

    can_scale = True
    w = 0.8
    default_aspect = w
    node_pinnames = ('s+', 's-', 'p+', 'p-')
    pins = {'s+': ('rx', 0.5 * w, 0.5),
            's-': ('rx', 0.5 * w, -0.5),
            'p+': ('lx', -0.5 * w, 0.5),
            'p-': ('lx', -0.5 * w, -0.5)}
    misc = {'pdot+': ('', 0.1 - 0.5 * w, 0.34),
            'sdot+': ('', 0.5 * w - 0.1, 0.34),
            'pdot-': ('', 0.1 - 0.5 * w, -0.34),
            'sdot-': ('', 0.5 * w - 0.1, -0.34),
            'link': ('', 0, 0.15),
            'plabel': ('lx', -0.7, 0),
            'slabel': ('rx', 0.7, 0),
            'label': ('', 0, 0.48)}

    def draw(self, link=True, **kwargs):

        if not self.check():
            return ''

        misc= self.tweak_pins(self.misc)

        centre = self.midpoint(self.nodes[0], self.nodes[3])
        label_pos = self.tf(centre, misc['label'][1:3])
        pdot_pos = self.tf(centre, misc['pdot+'][1:3])
        sdot_pos = self.tf(centre, misc['sdot+'][1:3])

        s = r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
            label_pos, 0.5, self.s, self.label(**kwargs))

        if link:
            # TODO: allow for rotation
            width = (sdot_pos - pdot_pos).x
            link_pos = self.tf(centre, misc['link'][1:3])

            s += r'  \draw[<->] ([shift=(45:%.2f)]%s) arc(45:135:%.2f);''\n' % (
                width / 2, link_pos, width / 2)

        core = ''
        if self.classname in ('TFcore', 'TFtapcore'):
            core = 'two'
            warn('Class ' + self.classname + ' is deprecated.  Use core attribute')

        if 'core' in self.opts:
            core = self.opts['core']

        if core in ('two', True):
            # Draw core with two lines
            q = self.tf(centre, ((-0.05, -0.2), (-0.05, 0.2),
                                 (0.05, -0.2), (0.05, 0.2)))
            s += self.draw_path(q[0:2], style='thick')
            s += self.draw_path(q[2:4], style='thick')
        elif core == 'one':
            # Draw core with one line
            q = self.tf(centre, ((0, -0.2), (0, 0.2)))
            s += self.draw_path(q[0:2], style='thick')
        elif core not in ('', False):
            raise ValueError('Invalid core attribute value ' + core)

        if self.turns:

            plabel_pos = self.tf(centre, misc['plabel'][1:3])
            slabel_pos = self.tf(centre, misc['slabel'][1:3])

            s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
                slabel_pos, 0.5, self.s, self.args[0])
            s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
                plabel_pos, 0.5, self.s, self.args[1])

        return s


class Transformer(TF1):
    """Transformer"""

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes[0:4]

        draw_args_str = self.draw_args_str(**kwargs)

        s = self.draw_cpt(n3.s, n4.s, cpt='inductor', args='scale=%s' %
                          self.scale, dargs=draw_args_str)
        s += self.draw_cpt(n2.s, n1.s, cpt='inductor', args='scale=%s' %
                          self.scale, dargs=draw_args_str)

        s += super(Transformer, self).draw(link=False, **kwargs)
        return s


class TFtap(TF1):
    """Transformer"""

    node_pinnames = ('s+', 's-', 'p+', 'p-', 'ptap', 'stap')
    w = 0.8
    pins = {'s+': ('rx', w, 1),
            's-': ('rx', w, 0),
            'p+': ('lx', 0, 1),
            'p-': ('lx', 0, 0),
            'ptap': ('lx', 0, 0.5),
            'stap': ('rx', w, 0.5)}

    @property
    def drawn_nodes(self):
        # Do not draw the taps.
        return self.nodes[0:4]

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes[0:4]

        draw_args_str = self.draw_args_str(**kwargs)

        s = self.draw_cpt(n3.s, n4.s, cpt='inductor', args='scale=%s' %
                          self.scale, dargs=draw_args_str)
        s += self.draw_cpt(n2.s, n1.s, cpt='inductor', args='scale=%s' %
                          self.scale, dargs=draw_args_str)

        s += super(TFtap, self).draw(link=False, **kwargs)
        return s


class TFp1s2(TF):
    """Transformer with 1 primary winding and 2 secondary windings"""

    can_invert = True
    can_mirror = True
    node_pinnames = ('s1+', 's1-', 's2+', 's2-', 'p1+', 'p1-')
    w = 0.8
    s = 0.4
    h = 1.0
    pins = {'s1+': ('rx', 0.5 * w, 0.5 * s + h),
             's1-': ('rx', 0.5 * w, 0.5 * s),
             's2+': ('rx', 0.5 * w, -0.5 * s),
             's2-': ('rx', 0.5 * w, -0.5 * s - h),
             'p1+': ('lx', -0.5 * w, 0.5 * h),
             'p1-': ('lx', -0.5 * w, -0.5 * h)}

    o = 0.25
    x = 0.5 * w + 0.15
    misc = {'sdot1+': ('', x, 0.5 * s + h - o),
            'sdot1-': ('', x, 0.5 * s + o),
            'sdot2+': ('', x, -0.5 * s - o),
            'sdot2-': ('', x, -0.5 * s - h + o),
            'pdot1+': ('', -x, 0.5 * h - o),
            'pdot1-': ('', -x, -0.5 * h + o),
            'plabel': ('lx', -x - o, 0),
            's2label': ('rx', x + o, 0.5 * (s + h)),
            's1label': ('rx', x + o, -0.5 * (s + h)),
            'label': ('', 0, 1.1)}

    @property
    def tpins(self):

        return self.tweak_pins(self.pins)

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4, n5, n6 = self.nodes[0:6]

        draw_args_str = self.draw_args_str(**kwargs)
        args = 'scale=%s' % self.scale
        if self.invert ^ self.mirror:
            args = ', '.join([args, 'mirror'])

        # n5 and n6 deliberately swapped for inductor symmetry
        s = self.draw_cpt(n5.s, n6.s, cpt='inductor', args=args,
                          dargs=draw_args_str)
        s += self.draw_cpt(n4.s, n3.s, cpt='inductor', args=args,
                           dargs=draw_args_str)
        s += self.draw_cpt(n2.s, n1.s, cpt='inductor', args=args,
                           dargs=draw_args_str)

        centre = (n2.pos + n3.pos + n5.pos + n6.pos) * 0.25
        if not self.nodots:

            misc= self.tweak_pins(self.misc)

            dots = {'TFptstt': '+++',
                    'TFptstb': '++-',
                    'TFptsbt': '+-+',
                    'TFptsbb': '+--'}[self.__class__.__name__]

            pdot1 = 'pdot1' + dots[0]
            sdot1 = 'sdot1' + dots[1]
            sdot2 = 'sdot2' + dots[2]
            pdot1_pos = self.tf(centre, misc[pdot1][1:3])
            sdot1_pos = self.tf(centre, misc[sdot1][1:3])
            sdot2_pos = self.tf(centre, misc[sdot2][1:3])

            s += r'  \draw (%s) node[circ] {};''\n' % pdot1_pos
            s += r'  \draw (%s) node[circ] {};''\n' % sdot1_pos
            s += r'  \draw (%s) node[circ] {};''\n' % sdot2_pos

        s += self.draw_core(centre, 0.9)
        label_pos = self.tf(centre, misc['label'][1:3])
        s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
            label_pos, 0.5, self.s, self.label(**kwargs))

        if self.turns:

            plabel_pos = self.tf(centre, misc['plabel'][1:3])
            s1label_pos = self.tf(centre, misc['s1label'][1:3])
            s2label_pos = self.tf(centre, misc['s2label'][1:3])

            s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
                s2label_pos, 0.5, self.s, self.args[0])
            s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
                s1label_pos, 0.5, self.s, self.args[1])
            s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
                plabel_pos, 0.5, self.s, self.args[2])
        return s

class TFp1s1(TF):
    """Transformer with 1 primary winding and 1 secondary winding"""

    can_invert = True
    can_mirror = True
    node_pinnames = ('s1+', 's1-', 'p1+', 'p1-')
    w = 0.8
    h = 1
    pins = {'s1+': ('rx', 0.5 * w, 0.5 * h),
             's1-': ('rx', 0.5 * w, -0.5 * h),
             'p1+': ('lx', -0.5 * w, 0.5 * h),
             'p1-': ('lx', -0.5 * w, -0.5 * h)}

    o = 0.15
    x = 0.5 * w + o
    y = 0.5 * h - o
    misc = {'pdot1+': ('', -x, y),
            'pdot1-': ('', -x, -y),
            'sdot1+': ('', x, y),
            'sdot1-': ('', x, -y),
            'plabel': ('lx', -x - o, 0),
            'slabel': ('rx', x + o, 0),
            'label': ('', 0, 0.48)}

    @property
    def tpins(self):

        return self.tweak_pins(self.pins)

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2, n3, n4 = self.nodes[0:4]

        draw_args_str = self.draw_args_str(**kwargs)
        args = 'scale=%s' % self.scale
        if self.invert ^ self.mirror:
            args = ', '.join([args, 'mirror'])

        # n3 and n4 deliberately swapped for inductor symmetry
        s = self.draw_cpt(n3.s, n4.s, cpt='inductor', args=args,
                          dargs=draw_args_str)
        s += self.draw_cpt(n2.s, n1.s, cpt='inductor', args=args,
                           dargs=draw_args_str)

        centre = (n1.pos + n2.pos + n3.pos + n4.pos) * 0.25
        if not self.nodots:

            misc= self.tweak_pins(self.misc)

            dots = {'TFptst': '++',
                    'TFptsb': '+-'}[self.__class__.__name__]

            pdot1 = 'pdot1' + dots[0]
            sdot1 = 'sdot1' + dots[1]
            pdot1_pos = self.tf(centre, misc[pdot1][1:3])
            sdot1_pos = self.tf(centre, misc[sdot1][1:3])

            s += r'  \draw (%s) node[circ] {};''\n' % pdot1_pos
            s += r'  \draw (%s) node[circ] {};''\n' % sdot1_pos

        s += self.draw_core(centre, 0.3)

        label_pos = self.tf(centre, misc['label'][1:3])
        s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
            label_pos, 0.5, self.s, self.label(**kwargs))

        if self.turns:

            plabel_pos = self.tf(centre, misc['plabel'][1:3])
            slabel_pos = self.tf(centre, misc['slabel'][1:3])

            s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
                slabel_pos, 0.5, self.s, self.args[0])
            s += r'  \draw (%s) node[minimum width=%.1f] (%s) {%s};''\n' % (
                plabel_pos, 0.5, self.s, self.args[1])

        return s


class K(TF1):
    """Mutual coupling"""

    def __init__(self, sch, namespace, name, cpt_type, cpt_id, string,
                 opts_string, node_names, keyword, *args):

        self.Lname1 = args[0]
        self.Lname2 = args[1]
        super(K, self).__init__(sch, namespace, name,
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
        nodes = [self.sch.nodes[n] for n in L2.explicit_node_names + L1.explicit_node_names]
        return nodes
