from numpy import dot
from .cpt import Cpt


class FixedCpt(Cpt):

    can_stretch = False

    @property
    def centre(self):
        # Look for centre pin.
        for node in self.nodes:
            if node.name.endswith('.mid'):
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

        return centre + dot((x * self.w, y * self.h), self.R(angle_offset)) * scale
