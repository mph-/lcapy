"""This module performs parsing of SPICE-like netlists.  It uses a
custom parser rather than lex/yacc to give better error messages.

Copyright 2015, 2016 Michael Hayes, UCECE

"""

import re

# Could use a script to generate parser and parsing tables if speed
# was important.

class Param(object):

    def __init__(self, name, base, comment):
        
        self.name = name
        self.base = base
        self.comment = comment
        self.baseclass = None

    def is_valid(self, string):

        if self.baseclass is None:
            return True
        return self.baseclass.is_valid(string)


class Rule(object):

    def __init__(self, cpt_type, classname, params, comment, pos):
        
        self.type = cpt_type
        self.classname = classname
        self.params = params
        self.comment = comment
        self.pos = pos

    def __repr__(self):

        return self.type + 'name ' + ' '.join(self.params)

    def syntax_error(self, error, string):

        raise ValueError('Syntax error: %s when parsing %s\nExpected format: %s' % (error, string, repr(self)))        

    def process(self, paramdir, string, fields, name, namespace):

        params = self.params
        if len(fields) > len(params):
            extra = ''
            if '(' in string:
                extra = ' (perhaps enclose expressions with parentheses in {})'
            self.syntax_error('Too many args' + extra, string)

        nodes = []
        args = []
        for m, param in enumerate(params):

            if m >= len(fields):
                # Optional argument
                if param[0] == '[':
                    break
                self.syntax_error('Missing arg %s' % param, string)

            if param[0] == '[':
                param = param[1:-1]

            field = fields[m]

            if paramdir[param].base in ('pin', 'node'):
                if field[0] == '.':
                    # Note name contains namespace
                    field = name + field
                else:
                    field = namespace + field
                nodes.append(field)
            elif paramdir[param].base != 'keyword':
                args.append(field)

        return tuple(nodes), args


class Parser(object):

    def __init__(self, cpts, grammar):
        """cpts is a module containing a class for each component
        grammar is a module defining the syntax of a netlist"""

        # A string defining the syntax for a netlist
        rules = grammar.rules
        # A string defining parameters
        params = grammar.params
        # A string defining delimiter characters
        delimiters = grammar.delimiters
        # A string defining comment characters
        self.comments = grammar.comments

        self.cpts = cpts
        self.paramdir = {}
        self.ruledir = {}
        
        for param in params.split('\n'):
            self._add_param(param)

        for rule in rules.split('\n'):
            self._add_rule(rule)

        cpts = sorted(self.ruledir.keys(), key=len, reverse=True)

        self.cpt_pattern = re.compile("(%s)([#_\w']+)?" % '|'.join(cpts))
        # strings in curly braces are expressions so do not split.
        # Fix if have {expr1} {expr2}
        self.param_pattern = re.compile('\{.*\}|".*"|[^%s]+' % delimiters)

    def _add_param(self, string):

        if string == '':
            return

        fields = string.split(':')
        paramname = fields[0]
        fields = fields[1].split(';')
        parambase = fields[0].strip()
        comment = fields[1].strip()
        
        self.paramdir[paramname] = Param(paramname, parambase, comment)

    def _add_rule(self, string):

        if string == '':
            return

        fields = string.split(':')
        cpt_classname = fields[0]
        fields = fields[1].split(';')
        string = fields[0].strip()
        comment = fields[1].strip()
        
        fields = string.split(' ')
        params = fields[1:]

        cpt_type = fields[0][0:-4]
        
        pos = None
        newparams = ()
        for m, param in enumerate(params):
            if param[0] == '[':
                param = param[1:-1]
            if param not in self.paramdir:
                raise ValueError('Unknown parameter %s for %s' % (param, string))
            if pos is None and self.paramdir[param].base == 'keyword':
                pos = m

        if cpt_type not in self.ruledir:
            self.ruledir[cpt_type] = ()
        self.ruledir[cpt_type] += (Rule(cpt_type, cpt_classname,
                                        params, comment, pos), )

    def _syntax_error(self, string):

        raise ValueError('%s\nExpected format: %s' % (repr(rule)))

    def parse(self, string, parent=None):
        """Parse string and create object"""

        string = string.strip()
        if string == '':
            return None

        if string[0] in self.comments:
            return None

        self.string = string
        fields = string.split(';', 1)

        fields = self.param_pattern.findall(fields[0])

        # Strip {}, perhaps should do with regexp.
        for m, field in enumerate(fields):
            if field[0] == '{':
                fields[m] = fields[m][1:-1]
        
        name = fields.pop(0)
        parts = name.split('.')
        namespace = ''
        if len(parts) > 1:
            namespace = '.'.join(parts[0:-1]) + '.'
        name = parts[-1]

        match = self.cpt_pattern.match(name)
        if match is None:
            raise ValueError('Unknown component for %s' % name)

        groups = match.groups()
        cpt_type, cpt_id = groups[0], groups[1]
        if cpt_id is None:
            cpt_id = ''

        # This is the most hackery aspect of this parser where we
        # choose the rule pattern based on a keyword.  If the
        # keyword is not present, default to first rule pattern.
        # Perhaps a factory should sort this out?
        rule = self.ruledir[cpt_type][0]
        keyword = ''
        for rule1 in self.ruledir[cpt_type]:
            pos = rule1.pos
            if pos is None:
                continue
            if len(fields) > pos and fields[pos].lower() == rule1.params[pos]:
                rule = rule1
                keyword = rule1.params[pos]
                break

        name = namespace + cpt_type + cpt_id
        if cpt_id == '' and parent is not None:
            name += parent._make_anon(cpt_type)

        nodes, args = rule.process(self.paramdir, string, fields, name, 
                                   namespace)

        fields = string.split(';')
        opts_string = fields[1].strip() if len(fields) > 1 else '' 

        keyword = (pos, keyword)
            
        return self.cpts.make(rule.classname, parent, name,
                              cpt_type, cpt_id, string, opts_string,
                              tuple(nodes), keyword, *args)
