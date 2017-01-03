import lcapy.grammar as grammar
from lcapy.parser import Parser


class NetfileMixin(object):

    def _init_parser(self, cpts):
        self.parser = Parser(cpts, grammar)        
        self.namespace = ''
        self._anon = {}

    def _make_anon(self, kind):
        """Make identifier for anonymous component"""

        if kind not in self._anon:
            self._anon[kind] = 0
        self._anon[kind] += 1        
        return 'anon' + str(self._anon[kind])

    def _include(self, string):

        parts = string.split(' ')
        if len(parts) < 2 or parts[0] != '.include':
            raise ValueError('Expecting include filename in %s' % string)
        filename = parts[1]
        if len(parts) == 2:
            return self._netfile_add(filename, self.namespace)
        
        if len(parts) != 4 and parts[2] != 'as':
            raise ValueError('Expecting include filename as name in %s' % string)
        name = parts[3]
        namespace = self.namespace
        self.namespace = name + '.' + namespace
        ret = self._netfile_add(filename, self.namespace)        
        self.namespace = namespace
        return ret

    def _parse(self, string, namespace=''):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        if string == '':
            return None            

        if string[0] == ';':
            if hasattr(self, 'opts'):
                self.opts.add(string[1:])
            return None

        if string[0:9] == '.include ':
            self._include(string)
            return None

        cpt = self.parser.parse(namespace + string, self)
        return cpt

    def add(self, string):
        """Add a component to the netlist.
        The general form is: 'Name Np Nm args'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        self._add(string)
        self._invalidate()

    def _add(self, string, namespace=''):
        """The general form is: 'Name Np Nm symbol'
        where Np is the positive node and Nm is the negative node.

        A positive current is defined to flow from the positive node
        to the negative node.
        """

        if '\n' in string:
            lines = string.split('\n')
            for line in lines:
                self._add(line.strip(), namespace)
            return

        cpt = self._parse(string, namespace)
        if cpt is not None:
            self._cpt_add(cpt)

    def _netfile_add(self, filename, namespace=''):
        """Add the nets from file with specified filename"""

        file = open(filename, 'r')

        lines = file.readlines()

        for line in lines:
            self._add(line, namespace)
