"""This module provides the SubNetlistSimplifyMixin class.

Copyright 2022--2023 Michael Hayes, UCECE

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
        subset_list.sort()

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

        newname = self.namer(name[0] + 'sim', self.elements)
        net1 = elt._new_value(total, ic)
        parts = net1.split(' ', 1)
        net1 = newname + ' ' + parts[1]

        net.add(net1)
        net.remove(name)

        for name1 in subset_list[1:]:
            # Replace with wire or open-circuit.
            if series:
                net1 = self.elements[name1]._netmake_W()
            else:
                net1 = self.elements[name1]._netmake_O()

            # Avoid creating open-circuit components.
            if True and series:
                net.add(net1)

            net.remove(name1)

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

        # TODO, remove dangling wires connected to the removed components.

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

        Note, if there are no circuits, e.g., a series network, then
        all the components will be removed.

        `select` is a list of component names to consider for simplification.
        If `None`, all components are considered.

        If `explain` is True, the reason for a simplification is printed.
        If `modify` is False, no modifications are performed.

        """

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

    def get_next_simplify_elements(self, series: bool = False, parallel: bool = False, debug: bool = False) -> list[str]:
        """
        The function returns two elements from a list. The list of elements is sorted by name so with the same
        components in a circuit the order of the returned elements will alway be the same. The elements in the returned
        list are either in series or in parallel depending on `series` and `parallel`.
        :param series: if true the function returns the first two elements which are in series.
        :param parallel: if true the function returns the first two elements which are in parallel.
        :param debug: prints the returned elements into a textfile named debug.txt
        :return: a list with two string elements
        """

        if series and parallel:
            raise AssertionError('You cannot specify both series and parallel. series || parallel = True')
        if not series and not parallel:
            raise AssertionError('You have to specify series or parallel. series || parallel = True')

        net = self.copy()

        if series and net.in_series():
            elements = list(net.in_series()[0])
        elif parallel and net.in_parallel():
            elements = list(net.in_parallel()[0])
        else:
            return []

        elements.sort()

        from lcapy import Circuit
        assert isinstance(net, Circuit)

        for element in elements:
            if element in net.components.voltage_sources \
                    or element in net.components.current_sources:
                elements.remove(element)

        if debug:
            print(elements[0:2])
            f = open("debug.txt", "a")
            f.write(str(elements[0:2])+"\n")

        return elements[0:2]

    def find_new_cpt_name(self, oldCpts: set, newCpts: set) -> str:
        """
        takes the components (cpts) of a Circuit before and after simplification and returns the name of the
        simplified element.
        :param oldCpts: set(Circuit.cpts) before Circuit.simplify()
        :param newCpts: set(Circuit.cpts) after Circuit.simplify()
        :return: string of the new entry
        """
        diffComp = newCpts - oldCpts
        diffComp = [component for component in diffComp if not component[0] == "W"]

        if diffComp:
            return diffComp[0]
        else:
            warn("No simplification performed, might not work as expected", RuntimeWarning)
            return ""

    def simplify_two_cpts(self, net, selected):
        """
        simplifies two selected componentes of a circuit and returns the new Circuit and the name of the simplified
        component.
        :param net: Circuit to be simplified
        :param selected: the thwo components to be simplified
        :return: the simplified circuit as a circuit and the name of the simplified component as a string
        """
        if len(selected) > 2:
            warn(f"first two components selected, length exceeded 2", RuntimeWarning)
            selected = selected[0:2]

        oldCpts = set(net.cpts)
        net = net.simplify(select=selected)
        newCpts = set(net.cpts)
        newCptName = self.find_new_cpt_name(oldCpts, newCpts)
        return net, newCptName

    def simplify_stepwise(self, limit: int = 100, debug: bool = False) -> list[tuple]:
        """
        Simplifies the circuit it is called on stepwise and returns a list of tupels which represent all steps
        that where made to simplify the circuit. The tuple contains the under return specified elements.
        Return value can be used as Input for Solution()-Object.
        :param limit: How many iterations are take before aborting the simplification process.
        :param debug: print debug info in debug.txt
        :return: a list of tupels with [(Circuit, StepComponent1, StepComponent2, StepComponentsCombined, Relation),...]
        """

        net = self.copy()
        steps = [(net, None, None, None, None)]

        if debug:
            f = open("debug.txt", "a")
            from datetime import datetime
            f.write(f"--- [{datetime.now()}] ---------\n")
            f.close()

        for i in range(0, limit):

            selected = net.get_next_simplify_elements(series=True, debug=debug)
            if len(selected) > 1:
                net, newCptName = self.simplify_two_cpts(net, selected=selected)
                steps.append((net, selected[0], selected[1], newCptName, "series"))
                continue

            selected = net.get_next_simplify_elements(parallel=True, debug=debug)
            if len(selected) > 1:
                net, newCptName = self.simplify_two_cpts(net, selected=selected)
                steps.append((net, selected[0], selected[1], newCptName, "parallel"))
                continue

            if net.in_series():
                lenSeries = len(list(net.in_series()[0]))
            else:
                lenSeries = 0

            if net.in_parallel():
                inParallel = list(net.in_parallel()[0])
                exclude = net.components.voltage_sources
                exclude.extend(net.components.current_sources)
                inParallel = [elem for elem in inParallel if elem not in exclude]
                lenParallel = len(inParallel)
            else:
                lenParallel = 0

            if lenSeries <= 1 and lenParallel <= 1:
                break

        return steps





    def simplify(self, select=None, ignore=None, passes=0, series=True,
                 parallel=True, dangling=False, disconnected=False,
                 explain=False, modify=True, keep_nodes=None):
        """Simplify a circuit by combining components in series and combining
        components in parallel.

        If `dangling` is True, then dangling components are removed.
        Dangling components are not in a circuit.

        If `disconnected` is True, then disconnected components are removed.
        Disconnected components are not connected to any other components.

        This performs a number of passes specified by `passes`.  If zero,
        this iterates until no more simplifications can be performed.

        `select` is a list of component names to consider for simplification.
        If `None`, all components are considered.

        If `explain` is True, the reason for a simplification is printed.
        If `modify` is False, no modifications are performed.

        See also `simplify_series`, `simplify_parallel`, `remove_dangling`,
        and `remove_dangling_wires`.

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
