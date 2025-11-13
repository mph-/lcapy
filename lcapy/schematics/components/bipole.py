from .stretchycpt import StretchyCpt
from ..utils import check_boolean
from ...label import Label
from ...latex import latex_format_label
from ...labelmaker import LabelMaker


class Bipole(StretchyCpt):
    """Bipole"""

    can_mirror = True
    can_invert = True
    can_scale = True

    node_pinnames = ('+', '-')
    aliases = {'p': '+', 'n': '-'}
    pins = {'+': ('lx', -0.5, 0),
            '-': ('rx', 0.5, 0)}
    misc = {'dot+': (-0.3, -0.1),
            'dot-': (0.3, -0.1)}
    # FIXME
    # auxiliary = {'mid': ('c', 0.0, 0.0)}

    def label_make(self, **kwargs):

        label = None
        annotation = None

        label_values = check_boolean(kwargs.get('label_values', True))
        label_ids = check_boolean(kwargs.get('label_ids', True))
        annotate_values = check_boolean(kwargs.get('annotate_values', False))
        style = kwargs.get('label_style', '=')

        if style == 'none':
            return label, annotation
        elif style == 'stacked':
            delimiter = '\\\\'
        elif style == 'aligned':
            delimiter = '='
        elif style == 'equals':
            delimiter = '='
        elif style == 'split':
            pass
        elif style == 'name':
            label_ids = True
            label_values = False
        elif style == 'value':
            label_ids = False
            label_values = True
        else:
            raise ValueError('Unknown label_style ' + style)

        label_value_style = kwargs.get('label_value_style', 'eng3')

        id_label, value_label = LabelMaker().make(self, style=label_value_style)

        id_label = latex_format_label(id_label)
        value_label = latex_format_label(value_label)

        # Generate default label.
        if (label_ids and label_values and id_label != ''
                and value_label != '' and id_label != value_label):
            if annotate_values or style == 'split':
                label = Label('l', id_label)
                annotation = Label('a', value_label)
            else:
                val = id_label + delimiter + value_label
                if '=' in val:
                    val = '{' + val + '}'
                label = Label('l', val)
        elif label_ids and id_label != '':
            label = Label('l', id_label)
        elif label_values and value_label != '':
            label = Label('l', value_label)

        return label, annotation

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n1, n2 = self.nodes[0:2]

        if self.wire:
            # With this option, draw component as a piece of wire.
            # This is useful for hiding the control voltage source
            # required for a CCVS or a CCCS.
            dargs = self.draw_args(self.opts, **kwargs)
            s = self.draw_wire(n1.s, n2.s, dargs)
            return s

        tikz_cpt = self.tikz_cpt
        if self.variable:
            if self.type in ('C', 'R', 'L'):
                tikz_cpt = 'v' + tikz_cpt
            else:
                raise ValueError('Component %s not variable' % self.name)
        if self.type == 'C' and self.args[1] is not None and 'v' not in self.opts:
            # Have initial value
            self.opts['v'] = 'v0{=}%s' % self.args[1]
        elif self.type == 'L' and self.args[1] is not None and 'f' not in self.opts:
            # Have initial value
            self.opts['f'] = 'i0{=}%s' % self.args[1]
        elif self.type in ('SW', 'SWno', 'SWnc') and 'l' not in self.opts:
            self.opts['l'] = 't{=}%s' % self.args[0]
        elif self.type == 'MISC':
            if self.kind is None:
                raise ValueError('Kind must be specified for %s' % self)
            tikz_cpt = self.kind
        else:
            if self.kind is not None:
                if self.kind not in self.kinds:
                    raise ValueError('Unknown kind %s for %s: known kinds %s'
                                     % (self.kind, self.name,
                                        ', '.join(self.kinds.keys())))
                tikz_cpt = self.kinds[self.kind]

            if self.style is not None:
                if self.style not in self.styles:
                    raise ValueError('Unknown style %s for %s: known styles %s'
                                     % (self.style, self.name,
                                        ', '.join(self.styles.keys())))
                tikz_cpt += self.styles[self.style]

        flip = check_boolean(kwargs.get('label_flip', False))

        if self.type in ('V', 'I', 'E', 'F', 'G', 'H', 'BAT'):
            # The node order has changed with different versions of
            # Circuitikz.  First there was the `old` order for
            # versions before 0.5.  Versions from 0.5 up to 0.9 use
            # what is called the `noold` order.  Versions from 0.9
            # introduce `RP` and `EP` orders.  These make more sense
            # but switch the node order for sources compared to the
            # `noold` order back to the `old` order.

            n1, n2 = n2, n1
            if self.left or self.down:
                # Draw label on LHS for vertical cpt and below
                # for horizontal cpt.
                flip = not flip

        cargs = self.cpt_args(self.opts, **kwargs)
        dargs = self.draw_args(self.opts, **kwargs)

        # Create default label if not overridden
        if not self.labels.label:
            label, annotation = self.label_make(**kwargs)
            if label:
                self.labels.label = label
            if annotation:
                self.labels.annotation = annotation

        # Add the specified labels
        cargs.extend(self.labels.args(flip))

        if self.mirror:
            cargs.append('mirror')
        if self.invert:
            cargs.append('invert')

        if self.scale != 1.0:
            cargs.append('bipoles/length=%.2fcm' % (
                self.sch.cpt_size * self.scale))

        cargs.append('n=' + self.s)

        s = self.draw_cpt(n1.s, n2.s, tikz_cpt, cargs, dargs)

        if self.type == 'L' and 'dot' in self.opts:
            dot = self.opts['dot']
            if dot != '':
                if dot not in ('+', '-'):
                    raise ValueError('Invalid value for dot attribute: ' + dot)
                centre = self.midpoint(self.nodes[0], self.nodes[1])
                dot_pos = self.tf(centre, self.misc['dot' + dot])
                s += r'  \draw (%s) node[circ] {};''\n' % dot_pos

        return s
