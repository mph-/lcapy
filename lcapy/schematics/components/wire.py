from warnings import warn
from .bipole import Bipole
from ...schemmisc import Steps


class Wire(Bipole):

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

        dargs = self.draw_args(self.opts, **kwargs)

        if self.steps is None:
            s = self.draw_wire(n1.s, n2.s, style=style,
                               startarrow=startarrow,
                               endarrow=endarrow,
                               dargs=dargs)
        else:
            steps = Steps(self.steps, n1.pos, n2.pos)

            s = self.draw_stepped_wire(n1.s, steps,
                                       style=style,
                                       startarrow=startarrow,
                                       endarrow=endarrow,
                                       dargs=dargs)

        if self.labels.voltage:
            # Well there can be an EMF if a changing magnetic flux passes
            # through the loop.
            warn('Ignoring voltage label; there is no voltage drop across an ideal wire!')

        return s
