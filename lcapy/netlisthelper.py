"""This module provides the NetlistHelper class used by the
NetlistMaker and LadderMaker classes.

Copyright 2020--2022 Michael Hayes, UCECE

"""

from .componentnamer import ComponentNamer


class NetlistHelper(object):

    evalf = False

    @property
    def _node(self):

        if not hasattr(self, '_node_counter'):
            self._node_counter = 0
        ret = self._node_counter
        self._node_counter += 1
        return ret

    def _make_nodes(self, *nodes):

        return [self._node if node is None else node for node in nodes]

    def _make_name(self, cpt_type, args):
        """Make component name"""

        if not hasattr(self, '_namer'):
            self._namer = ComponentNamer()

        name = None
        # See if argument is suitable for the cpt_id
        if (len(args) > 0 and isinstance(args[0], str) and
                args[0].startswith(cpt_type) and args[0].isalnum()):
            name = args[0]

            # If have already used this name, then resort to auto name.
            if name in self.names:
                name = None

        if name is None:
            name = self._namer.name(cpt_type)

        self.names.append(name)

        return name

    def _netarg(self, arg):

        if self.evalf:
            try:
                # The arg may not be an Expr.
                arg = arg.evalf(n=self.evalf)
            except:
                pass

        try:
            # This is primarily for Superposition values.
            argstr = arg._netrepr()
        except:
            argstr = str(arg)

        # TODO: make more robust to catch expressions.
        if (('(' in argstr) or (')' in argstr) or (' ' in argstr) or
                (',' in argstr) or ('*' in argstr) or ('/' in argstr)):
            return '{%s}' % argstr
        return argstr

    def _netargs(self, net):

        return ' '.join([self._netarg(arg) for arg in net.args])
