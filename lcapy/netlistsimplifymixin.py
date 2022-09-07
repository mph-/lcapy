"""This module provides the SubNetlistSimplifyMixin class.

Copyright 2022 Michael Hayes, UCECE

"""

from .expr import expr
from warnings import warn
from collections import OrderedDict


class NetlistSimplifyMixin:

    def _do_simplify_combine(self, string, subset, net,
                             explain=False, add=False, series=False):

        if explain:
            print(string % subset)

        subset_list = list(subset)

        if add:
            total = expr(0)
            for name in subset_list:
                total += expr(self.elements[name].cpt.args[0])
        else:
            total = expr(0)
            for name in subset_list:
                total += (1 / expr(self.elements[name].cpt.args[0]))
            total = 1 / total

        if explain:
            print('%s combined value = %s' % (subset, total))

        ic = None
        name = subset_list[0]
        elt = self.elements[name]
        if elt.cpt.has_ic:
            ic = expr(0)
            for name1 in subset_list:
                ic += expr(self.elements[name1].cpt.args[1])

            if explain:
                print('%s combined IC = %s' % (subset, ic))

        newname = self.namer(name[0] + 't', self.elements)
        net1 = elt._new_value(total, ic)
        parts = net1.split(' ', 1)
        net1 = newname + ' ' + parts[1]

        elt = self._parse(net1)

        # Overwrite component with one having total value.  _fixup will
        # fix element keys later on, once all simplifications are performed.
        net.elements[name] = elt

        for name1 in subset_list[1:]:
            # Replace with wire or open-circuit.
            if series:
                net1 = self.elements[name1]._netmake_W()
            else:
                net1 = self.elements[name1]._netmake_O()
            # Remove component
            elt = self._parse(net1)
            net.elements[name1] = elt

        return True

    def _check_ic(self, subset):

        subset = subset.copy()
        name = subset.pop()
        has_ic = self.elements[name].has_ic

        okay = True
        for name1 in subset:
            if self.elements[name1].has_ic != has_ic:
                warn('Incompatible initial conditions for %s and %s' %
                     (name, name1))
                okay = False
        if not has_ic:
            return okay
        ic = self.elements[name].cpt.args[1]
        for name1 in subset:
            if self.elements[name1].cpt.args[1] != ic:
                warn('Incompatible initial conditions for %s and %s' %
                     (name, name1))
                okay = False

        return okay

    def _fixup(self):
        """Rename keys to fix things up for removed components."""

        newelements = OrderedDict()
        for k, v in self.elements.items():
            newelements[v.name] = v
        self._elements = newelements

    def _simplify_combine_series(self, skip, explain=False):

        net = self.copy()
        changed = False

        for aset in net.in_series():
            aset -= skip
            subsets = net._find_combine_subsets(aset)
            for k, subset in subsets.items():
                if k == 'I':
                    warn('Netlist has current sources in series: %s' % subset)
                elif k in ('R', 'NR', 'L', 'V', 'Z'):
                    if k == 'L' and not self._check_ic(subset):
                        continue
                    changed |= self._do_simplify_combine('Can add in series: %s',
                                                         subset, net, explain, True, True)
                elif k in ('C', 'Y'):
                    changed |= self._do_simplify_combine('Can combine in series: %s',
                                                         subset, net, explain, False, True)
                else:
                    raise RuntimeError('Internal error')

        if changed:
            net._fixup()

        return net, changed

    def _simplify_combine_parallel(self, skip, explain=False):

        net = self.copy()
        changed = False

        for aset in net.in_parallel():
            aset -= skip
            subsets = net._find_combine_subsets(aset)
            for k, subset in subsets.items():
                if k == 'V':
                    warn('Netlist has voltage sources in parallel: %s' % subset)
                elif k in ('R', 'NR', 'L', 'Z'):
                    changed |= self._do_simplify_combine('Can combine in parallel: %s',
                                                         subset, net, explain, False, False)
                elif k in ('C', 'Y', 'I'):
                    if k == 'C' and not self._check_ic(subset):
                        continue
                    changed |= self._do_simplify_combine('Can add in parallel: %s',
                                                         subset, net, explain, True, False)
                else:
                    raise RuntimeError('Internal error')

        if changed:
            # TODO, remove dangling wires connected to the removed components.
            net._fixup()

        return net, changed

    def _simplify_redundant_series(self, skip, explain=False):

        net = self.copy()
        changed = False

        for aset in net.in_series():
            Iname = None
            for name in aset:
                cpt = self._elements[name]
                if cpt.type == 'I':
                    Iname = name
                    break
            if Iname is not None:
                for name in aset:
                    cpt = self._elements[name]
                    if cpt.type != 'I':
                        warn('Have redundant %s in series with %s' %
                             (name, Iname))

        return net, False

    def _simplify_redundant_parallel(self, skip, explain=False):

        net = self.copy()
        changed = False

        for aset in net.in_parallel():
            Vname = None
            for name in aset:
                cpt = self._elements[name]
                if cpt.type == 'V':
                    Vname = name
                    break
            if Vname is not None:
                for name in aset:
                    cpt = self._elements[name]
                    if cpt.type != 'V':
                        warn('Have redundant %s in parallel with %s' %
                             (name, Vname))

        return net, False

    def _keep_dangling(self, cpt, keep_nodes):

        for node in cpt.nodes:
            if node.is_dangling and node.name in keep_nodes:
                return True

        return False

    def _remove_dangling(self, skip, explain=False, keep_nodes=None):

        new = self._new()
        changed = False

        for cpt in self._elements.values():
            if (cpt.is_dangling and cpt.name not in skip
                    and not self._keep_dangling(cpt, keep_nodes)):
                if explain:
                    print('Removing dangling component %s' % cpt.name)
                changed = True
            else:
                new._add(cpt._copy())

        return new, changed

    def _remove_disconnected(self, skip, explain=False, keep_nodes=None):

        new = self._new()
        changed = False

        for cpt in self._elements.values():
            if (cpt.is_disconnected and cpt.name not in skip
                    and not self._keep_dangling(cpt, keep_nodes)):
                if explain:
                    print('Removing disconnected component %s' % cpt.name)
                changed = True
            else:
                new._add(cpt._copy())

        return new, changed

    def _simplify_series(self, skip, explain=False):

        net, changed = self._simplify_redundant_series(skip, explain)
        net, changed2 = net._simplify_combine_series(skip, explain)
        return net, changed or changed2

    def _simplify_parallel(self, skip, explain=False):

        net, changed = self._simplify_redundant_parallel(skip, explain)
        net, changed2 = net._simplify_combine_parallel(skip, explain)
        return net, changed or changed2

    def remove_dangling(self, select=None, ignore=None, passes=0, explain=False,
                        modify=True, keep_nodes=None):
        """Simplify a circuit by removing dangling components.

        This performs a number of passes specified by `passes`.  If zero,
        this iterates until no more simplifications can be performed.

        `select` is a list of component names to consider for simplification.
        If `None`, all components are considered.

        If `explain` is True, the reason for a simplification is printed.
        If `modify` is False, no modifications are performed."""

        return self.simplify(select=select, ignore=ignore, passes=passes,
                             explain=explain, modify=modify,
                             series=False, parallel=False, dangling=True,
                             keep_nodes=keep_nodes)

    def remove_dangling_wires(self, passes=0, explain=False,
                              modify=True, keep_nodes=None):
        """Simplify a circuit by removing dangling wires.

        This also removes disconnected open-circuit components
        after dangling wires have been removed.

        This performs a number of passes specified by `passes`.  If zero,
        this iterates until no more simplifications can be performed.

        If `explain` is True, the reason for a simplification is printed.
        If `modify` is False, no modifications are performed."""

        new = self.simplify(select=self.components.wires,
                            passes=passes, explain=explain,
                            modify=modify, series=False,
                            parallel=False, dangling=True, disconnected=False,
                            keep_nodes=keep_nodes)

        # Remove disconnected open-circuit components.
        return new.simplify(select=self.components.open_circuits,
                            passes=passes, explain=explain,
                            modify=modify, series=False,
                            parallel=False, dangling=False, disconnected=True,
                            keep_nodes=keep_nodes)

    def remove_disconnected(self, select=None, ignore=None, passes=0,
                            explain=False, modify=True, keep_nodes=None):
        """Simplify a circuit by removing disconnected components.

        This performs a number of passes specified by `passes`.  If zero,
        this iterates until no more simplifications can be performed.

        `select` is a list of component names to consider for simplification.
        If `None`, all components are considered.

        If `explain` is True, the reason for a simplification is printed.
        If `modify` is False, no modifications are performed."""

        return self.simplify(select=select, ignore=ignore, passes=passes,
                             explain=explain, modify=modify,
                             series=False, parallel=False, dangling=False,
                             disconnected=True, keep_nodes=keep_nodes)

    def simplify_series(self, select=None, ignore=None, passes=0,
                        explain=False, modify=True, keep_nodes=None):
        """Simplify a circuit by combining components in series.

        This performs a number of passes specified by `passes`.  If zero,
        this iterates until no more simplifications can be performed.

        `select` is a list of component names to consider for simplification.
        If `None`, all components are considered.

        If `None`, all components are considered.
        If `explain` is True, the reason for a simplification is printed.
        If `modify` is False, no modifications are performed."""

        return self.simplify(select=select, ignore=ignore, passes=passes,
                             explain=explain, modify=modify,
                             series=True, parallel=False, dangling=False,
                             keep_nodes=keep_nodes)

    def simplify_parallel(self, select=None, ignore=None, passes=0,
                          explain=False, modify=True, keep_nodes=None):
        """Simplify a circuit by combining components in parallel.

        This performs a number of passes specified by `passes`.  If zero,
        this iterates until no more simplifications can be performed.

        `select` is a list of component names to consider for simplification.
        If `None`, all components are considered.

        If `None`, all components are considered.
        If `explain` is True, the reason for a simplification is printed.
        If `modify` is False, no modifications are performed."""

        return self.simplify(select=select, ignore=ignore, passes=passes,
                             explain=explain, modify=modify,
                             series=False, parallel=True, dangling=False,
                             keep_nodes=keep_nodes)

    def simplify(self, select=None, ignore=None, passes=0, series=True,
                 parallel=True, dangling=False, disconnected=True,
                 explain=False, modify=True, keep_nodes=None):
        """Simplify a circuit by combining components in series, combining
        components in parallel, and removing dangling components.

        This performs a number of passes specified by `passes`.  If zero,
        this iterates until no more simplifications can be performed.

        `select` is a list of component names to consider for simplification.
        If `None`, all components are considered.

        If `explain` is True, the reason for a simplification is printed.
        If `modify` is False, no modifications are performed.

        See also `simplify_series`, `simplify_parallel`, `simplify_dangling`,
        and `simplify_dangling_wires`.
        """

        if keep_nodes is None:
            if '0' in self.nodes:
                keep_nodes = ['0']
            else:
                keep_nodes = []

        keep_nodes = [str(node) for node in keep_nodes]

        skip = set()

        if select is not None:
            skip = set(self._elements) - set(select)

        if ignore is not None:
            skip = skip.union(set(ignore))

        # Perhaps use num cpts?
        if passes == 0:
            passes = 100

        net = self
        for m in range(passes):
            changed = False
            if dangling:
                net, changed1 = net._remove_dangling(skip, explain, keep_nodes)
                changed = changed or changed1
            if disconnected:
                net, changed1 = net._remove_disconnected(
                    skip, explain, keep_nodes)
                changed = changed or changed1
            if series:
                net, changed1 = net._simplify_series(skip, explain)
                changed = changed or changed1
            if parallel:
                net, changed1 = net._simplify_parallel(skip, explain)
                changed = changed or changed1

            if not changed:
                break
        if not modify:
            return self

        return net
