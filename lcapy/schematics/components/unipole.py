from .cpt import Cpt


class Unipole(Cpt):

    place = False

    node_pinnames = ('+', )
    aliases = {'p': '+'}
    pins = {'+': ('lx', 0, 0)}

    def draw(self, **kwargs):

        if not self.check():
            return ''

        n = self.nodes[0]
        q = self.tf(n.pos, ((self.xoffset, self.yoffset)))
        tikz_cpt = self.tikz_cpt

        if self.kind is not None:
            if self.kind not in self.kinds:
                raise ValueError('Unknown kind %s for %s: known kinds %s'
                                 % (self.kind, self.name,
                                    ', '.join(self.kinds.keys())))
            tikz_cpt = self.kinds[self.kind]

        xscale = self.scale
        yscale = self.scale
        if self.mirror:
            yscale = -yscale
        if self.invert:
            xscale = -xscale

        args = self.draw_args(self.opts, **kwargs)
        if self.angle != 0:
            args.append('rotate=%d' % self.angle)
        if xscale != 1:
            args.append('xscale=%f' % xscale)
        if yscale != 1:
            args.append('yscale=%f' % yscale)
        anchor = self.anchor_opt(self)
        if anchor is not None:
            args.append('anchor=' + anchor)

        label = self.label(**kwargs)
        label = self.label_tweak(label, xscale, yscale, self.angle)
        s = self.draw_cptnode(q, tikz_cpt, args, '', label)
        return s
