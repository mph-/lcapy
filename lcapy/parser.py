"""This module performs parsing of SPICE-like netlists.  It uses a
custom parser rather than lex/yacc to give better error messages.

Copyright 2015--2022 Michael Hayes, UCECE

"""

import re

# Could use a script to generate parser and parsing tables if speed
# was important.

# Each line of a netlist (called net here) has a name followed by a
# number of required parameters (params), a number of optional
# parameters, an optional semicolon, and optional comma separated
# key-value pairs.
#
# Each param can be a node, keyword, name, or value.
#
# Here's an example: I1 4 5 {Piecewise((V/L, t >= 0))}; right
# The name is I1; 4, 5 are nodes; the expression is brackets is a value.
#
# The name or value params are referred to as args.


def split(s, delimiters):
    """Split string by specified delimiters but not if a delimiter is
    within curly brackets {} or ""."""

    parts = []
    current = []
    close_bracket = ''
    bracket_stack = []
    for c in (s + delimiters[0]):
        if c in delimiters and len(bracket_stack) == 0:
            if len(current) > 0:
                parts.append(''.join(current))
            current = []
        else:
            if c == close_bracket:
                close_bracket = bracket_stack.pop()
            elif c == '{':
                bracket_stack.append(close_bracket)
                close_bracket = '}'
            elif c == '"':
                bracket_stack.append(close_bracket)
                close_bracket = '"'
            current.append(c)
    if close_bracket != '':
        raise ValueError('Missing %s in %s' % (close_bracket, s))
    return parts


class ParamDef:

    def __init__(self, name, base, comment):

        self.name = name
        self.base = base
        self.comment = comment
        self.baseclass = None

    def is_valid(self, string):

        if self.baseclass is None:
            return True
        return self.baseclass.is_valid(string)


class Param:

    def __init__(self, paramstr, paramdict):

        self.optional = paramstr[0] == '['

        if self.optional:
            paramstr = paramstr[1:-1]

        parts = paramstr.split('=')
        self.name = parts[0]

        if self.name not in paramdict:
            raise ValueError('Unknown parameter ' + self.name)
        self.kind = paramdict[self.name].base

        if len(parts) == 1:
            self.default = None
        else:
            self.default = parts[1]

        self.lowercase_name = self.name.lower()


class Value:

    def __init__(self, param, name):
        """Set default value."""

        self.name = param.name
        self.assigned = False

        if param.default == 'name':
            default = name
        else:
            default = param.default
        self.default = default
        self.value = default


class Arg:

    def __init__(self, arg):

        if arg[0] in '{"':
            arg = arg[1:-1]
        self.value = arg


class Args(dict):
    pass


class Rule:

    def __init__(self, cpt_type, classname, fields, params, comment, pos):

        self.type = cpt_type
        self.classname = classname
        self.fields = fields
        self.params = params
        self.comment = comment
        self.pos = pos

    def __repr__(self):

        return self.type + 'name ' + ' '.join(self.fields[1:])

    def syntax_error(self, error, string):

        raise ValueError('Syntax error: %s when parsing %s\nExpected format: %s' % (
            error, string, repr(self)))

    def extract_nodes(self, string, fields, name, namespace):

        ruleparams = self.params
        nodes = []

        for m, ruleparam in enumerate(ruleparams):

            if ruleparam.kind not in ('pin', 'node'):
                continue

            if m >= len(fields):
                self.syntax_error('Missing node %s' % ruleparam.name, string)

            field = fields[m]
            if field[0] == '.':
                # Note, name contains namespace
                field = name + field
            else:
                field = namespace + field
            nodes.append(field)

        return tuple(nodes)

    def process(self, string, fields, name, namespace):

        ruleparams = self.params
        if len(fields) > len(ruleparams):
            extra = ''
            if '(' in string:
                extra = ' (perhaps enclose expressions with parentheses in {})'
            self.syntax_error('Too many args' + extra, string)

        nodes = self.extract_nodes(string, fields, name, namespace)

        args = []
        for m, ruleparam in enumerate(ruleparams):

            if ruleparam.kind not in ('name', 'value'):
                continue

            if m >= len(fields):
                if ruleparam.optional:
                    break
                self.syntax_error('Missing arg %s' % ruleparam, string)

            field = fields[m]
            args.append(field)

        return nodes, args

    def parse_args(self, net, name, args):

        args_dict = Args()

        values = []
        names = []

        # Determine parameters and default values
        for param in self.params:
            if (param.kind != 'value' and param.kind != 'name'):
                continue
            value = Value(param, name)

            values.append(value)
            names.append(value.name)
            args_dict[value.name] = value.default

        args_list = []

        # Handle unnamed params
        m = 0
        for arg in args:
            parts = split(arg, '=')
            if len(parts) > 1:
                break
            args_list.append(arg)
            values[m].assigned = True
            value = Arg(arg).value
            values[m].value = value
            args_dict[values[m].name] = value
            m += 1

        # Handle named params
        for arg in args[m:]:
            parts = split(arg, '=')
            if len(parts) < 2:
                self.syntax_error(
                    'Cannot have value %s after named param' % arg, net)
            parts = split(arg, '=')
            index = names.index(parts[0])
            if index < 0:
                self.syntax_error('Unknown param ' + parts[0], net)
            if values[index].assigned:
                self.syntax_error('Value %s already assigned' % parts[0], net)
            values[index].assigned = True
            arg = parts[1]
            args_dict[parts[0]] = Arg(arg).value

        args_list = []
        ignore = True
        for value in values[::-1]:
            if ignore and not value.assigned:
                continue
            ignore = False
            args_list.insert(0, value.value)

        return args_dict, args_list


class Parser:

    def __init__(self, cpts, grammar, allow_anon=False):
        """cpts is a module containing a class for each component
        grammar is a module defining the syntax of a netlist"""

        # A string defining the syntax for a netlist
        rules = grammar.rules
        # A string defining parameters
        params = grammar.params
        # A string defining delimiter characters
        self.delimiters = grammar.delimiters
        # A string defining comment characters
        self.comments = grammar.comments
        self.allow_anon = allow_anon

        self.cpts = cpts
        self.paramdict = {}
        self.ruledict = {}

        for param in params.split('\n'):
            self._add_param(param)

        for rule in rules.split('\n'):
            self._add_rule(rule)

        cpts = sorted(self.ruledict.keys(), key=len, reverse=True)

        # The symbol name must be a valid Sympy symbol name so
        # it cannot include symbols such as + and -.
        self.cpt_pattern = re.compile("(%s)([#_\w'?]+)?" % '|'.join(cpts))

    def _add_param(self, string):

        if string == '':
            return

        fields = string.split(':')
        paramname = fields[0]
        fields = fields[1].split(';', 1)
        parambase = fields[0].strip()
        comment = fields[1].strip()

        self.paramdict[paramname] = ParamDef(paramname, parambase, comment)

    def _add_rule(self, string):

        if string == '':
            return

        fields = string.split(':')
        cpt_classname = fields[0]
        fields = fields[1].split(';', 1)
        string = fields[0].strip()
        comment = fields[1].strip()

        fields = string.split(' ')

        # Skip the name part in the rule, e.g., only consider D from Dname.
        cpt_type = fields[0][0:-4]

        pos = None
        params = []
        for m, paramstr in enumerate(fields[1:]):

            param = Param(paramstr, self.paramdict)
            params.append(param)
            if pos is None and param.kind == 'keyword':
                pos = m

        if cpt_type not in self.ruledict:
            self.ruledict[cpt_type] = ()
        self.ruledict[cpt_type] += (Rule(cpt_type, cpt_classname, fields,
                                         params, comment, pos), )

    def parse(self, string, namespace='', parent=None):
        """Parse string and create object"""

        directive = False
        net = string.strip()
        if net == '':
            directive = True
        elif net[0] in self.comments:
            directive = True
        elif net[0] == ';':
            directive = True
        elif net[0] == '.':
            directive = True

        if directive:
            cpt_type = 'XX'
            cpt_id = ''
            name = 'XX'
            name += parent._make_anon_cpt_id(cpt_type)
            defname = namespace + cpt_type + cpt_id

            if string.startswith(';') and not string.startswith(';;'):
                opts_string = string[1:]
            else:
                opts_string = ''

            return self.cpts.make('XX', parent, '', defname, name,
                                  cpt_type, cpt_id, string, opts_string, (), '',
                                  Args())

        net = namespace + net
        parts = net.split(';', 1)

        fields = split(parts[0], self.delimiters)

        name = fields.pop(0)
        parts = name.split('.')
        namespace = ''
        if len(parts) > 1:
            namespace = '.'.join(parts[0:-1]) + '.'
        name = parts[-1]

        match = self.cpt_pattern.match(name)
        if match is None:
            raise ValueError(
                'Unknown component %s while parsing "%s"' % (name, net))

        groups = match.groups()
        cpt_type, cpt_id = groups[0], groups[1]
        if cpt_id is None:
            cpt_id = ''

        # This is the most hackery aspect of this parser where we
        # choose the rule pattern based on a keyword.  If the
        # keyword is not present, default to first rule pattern.
        # Perhaps a factory should sort this out?
        rule = self.ruledict[cpt_type][0]
        keyword = ''

        for rule1 in self.ruledict[cpt_type]:
            pos = rule1.pos
            if pos is None:
                continue
            if len(fields) > pos and fields[pos].lower() == rule1.params[pos].lowercase_name:
                rule = rule1
                keyword = rule1.params[pos]
                break

        defname = namespace + cpt_type + cpt_id
        name = defname
        if (cpt_id == '' and parent is not None
                and (cpt_type in ('A', 'W', 'O', 'P')) or self.allow_anon):
            name += parent._make_anon_cpt_id(cpt_type)
        elif cpt_id == '?':
            # Automatically name cpts to ensure they are unique
            name = name[:-1] + parent._make_anon_cpt_id(cpt_type)

        nodes, args = rule.process(net, fields, name, namespace)

        default_value = cpt_type + cpt_id
        args_dict, args_list = rule.parse_args(net, default_value, args)

        parts = net.split(';', 1)
        opts_string = parts[1].strip() if len(parts) > 1 else ''

        keyword = (pos, keyword)

        # self.cpts is either the mnacpts or schematic module
        return self.cpts.make(rule.classname, parent, namespace,
                              defname, name, cpt_type, cpt_id, net,
                              opts_string, tuple(nodes), keyword,
                              args_dict, *args_list)
