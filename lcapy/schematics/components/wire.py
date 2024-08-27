from warnings import warn
from .bipole import Bipole
from ...schemmisc import Steps


def arrow_map(name):

    try:
        name = {'tee': '|', 'otri': 'open triangle 60',
                'tri': 'triangle 60'}[name]
    except KeyError:
        pass
    return name


class Wire(Bipole):

    def draw_short(self, pos1, pos2, startarrow='', endarrow='',
                   style='', **kwargs):
        """Create a string to draw a circuitikz short between positions `pos1`
        and `pos2`.  `dargs` is a list or string of the draw options.

        The general form of the generated string is:

        \draw[dargs] (pos1) to [short, args] (pos2);

        Unlike a wire, a short can have labels.  However, it cannot
        have a startarrow.

        """

        if startarrow != '':
            warn('Ignoring start arrow.  If you want it, remove the label.')

        endarrow = arrow_map(endarrow)
        arrows = '-' + endarrow
        if arrows == '-':
            arrows = ''

        dargs = self.draw_args(self.opts, **kwargs)
        dargs.append(arrows)
        args = self.cpt_args(self.opts, **kwargs)
        args.extend(self.labels.args())

        s = self.draw_cpt(pos1, pos2, 'short', args, dargs)
        return s

    def draw_wire(self, pos1, pos2, startarrow='', endarrow='',
                  style='', **kwargs):
        """Create a string to draw a circuitikz wire between positions `pos1`
        and `pos2`.  `dargs` is a list or string of the draw options.

        The general form of the generated string is:

        \draw[-, dargs] (pos1) to (pos2);

        A wire cannot have labels; use a short instead.

        """

        startarrow = arrow_map(startarrow)
        endarrow = arrow_map(endarrow)

        arrows = startarrow + '-' + endarrow
        if arrows == '-':
            arrows = ''

        # Also use cpt args such as color
        dargs = self.draw_args(self.opts, **kwargs) + self.cpt_args(self.opts, **kwargs)
        dargs = [] if dargs is None else dargs
        dargs.insert(0, arrows)
        dargs.append(style)
        dargs = ', '.join([arg for arg in dargs if arg != ''])

        s = self.draw_cpt(pos1, pos2, '', dargs, '')
        return s

    def draw_stepped_wire(self, pos1, steps, startarrow='',
                          endarrow='', style='', **kwargs):

        path = '(%s)' % pos1
        for pos in steps:
            path += ' to (%s)' % pos

        startarrow = arrow_map(startarrow)
        endarrow = arrow_map(endarrow)
        arrows = startarrow + '-' + endarrow

        dargs = self.draw_args(self.opts, **kwargs)
        dargs = [] if dargs is None else dargs
        dargs.insert(0, arrows)
        dargs.append(style)
        dargs = ', '.join([arg for arg in dargs if arg != ''])

        s = r'  \draw[%s] %s;''\n' % (dargs, path)
        return s

    def draw(self, **kwargs):

        if not self.check():
            return ''

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

        if self.labels.voltage:
            # Well there can be an EMF if a changing magnetic flux passes
            # through the loop but if wire ideal will have an infinite current!
            warn('Ignoring voltage label; there is no voltage drop across an ideal wire!')

        if self.steps is None:

            label_args = self.labels.args()
            implicit = self.implicit_key(self.opts)

            # Don't draw labell if have implicit wire since
            # the label is used for the node.
            if label_args != [] and not implicit:
                s = self.draw_short(n1.s, n2.s, style=style,
                                    startarrow=startarrow,
                                    endarrow=endarrow, **kwargs)
            else:
                s = self.draw_wire(n1.s, n2.s, style=style,
                                   startarrow=startarrow,
                                   endarrow=endarrow, **kwargs)
        else:
            steps = Steps(self.steps, n1.pos, n2.pos)

            s = self.draw_stepped_wire(n1.s, steps, style=style,
                                       startarrow=startarrow,
                                       endarrow=endarrow, **kwargs)

        return s
