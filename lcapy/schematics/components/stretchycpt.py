from numpy import dot
from .cpt import Cpt


class StretchyCpt(Cpt):

    can_stretch = True

    def xtf(self, centre, offset, angle_offset=0.0):
        """Transform x coordinate."""

        # Note the size attribute is not used.   This only scales the x coords.
        if isinstance(offset[0], tuple):
            return [self.xtf(centre, offset1, angle_offset) for offset1 in offset]

        return centre + dot((offset[0] * self.w * self.scale, offset[1] * self.h), self.R(angle_offset)) * self.sch.node_spacing
