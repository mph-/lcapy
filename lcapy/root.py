from .quadraticroot import QuadraticRoot


class Root:

    def __init__(self, expr, n, damping=None):

        self.n = n

        d = None
        if damping is not None:
            # Perhaps parse as quadratic root, cubic root, etc.?
            d = QuadraticRoot.from_expr(expr, damping)

        if d is None:
            self.expr = expr
            self.conj = expr.conjugate()
        else:
            self.expr = d.expr
            self.conj = d.conjugate()

        self.damping = damping
        self._d = d

    def conjugate(self):
        return self.conj

    def is_conjugate_pair(self, root):

        return self.expr == root.conj


def pair_conjugates(roots_dict, damping=None):
    """Return dictionary of conjugate root pairs and a dictionary of the
    remaining single roots."""

    root_single_dict = roots_dict.copy()
    root_pair_dict = {}

    root_list = list(roots_dict)

    P = {}
    for root in root_list:
        P[root] = Root(root, 1, damping=damping)

    for i, root in enumerate(root_list):

        p = P[root]

        for root_c in root_list[i + 1:]:

            pc = P[root_c]
            if p.is_conjugate_pair(pc):

                root_single_dict.pop(root, None)
                root_single_dict.pop(root_c, None)

                o1 = roots_dict[root]
                o2 = roots_dict[root_c]
                if o1 == o2:
                    root_pair_dict[root, root_c] = o1
                elif o1 > o2:
                    root_pair_dict[root, root_c] = o2
                    root_single_dict[root] = o1 - o2
                else:
                    root_pair_dict[root, root_c] = o1
                    root_single_dict[root_c] = o2 - o1

    return root_pair_dict, root_single_dict
