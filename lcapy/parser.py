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

    def process(self, cpts, paramdir, name, string, fields):

        params = self.params
        if len(fields) > len(params):
            self.syntax_error('Too many args', string)

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
            if paramdir[param].base == 'node':
                if field.find('.') != -1:
                    self.syntax_error('Found . in node name %s' % field, string)
                nodes.append(field)
            elif paramdir[param].base != 'keyword':
                args.append(field)

        # Create instance of component object
        try:
            newclass = getattr(cpts, self.classname)
        except:
            newclass = cpts.classes[self.classname]

        cpt = newclass(name, string, tuple(nodes), *args)
        # Add named attributes for the args?   Lname1, etc.
        return cpt


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
        self.anon = {}
        self.paramdir = {}
        self.ruledir = {}
        
        for param in params.split('\n'):
            self._add_param(param)

        for rule in rules.split('\n'):
            self._add_rule(rule)

        cpts = self.ruledir.keys()
        cpts.sort(key=len, reverse=True)

        self.cpt_pattern = re.compile("(%s)([#_\w']+)?" % '|'.join(cpts))
        # strings in curly braces are expressions so do not split.
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

    def _anon_cpt_id(self, cpt_type):

        if cpt_type not in self.anon:
            self.anon[cpt_type] = 0

        self.anon[cpt_type] += 1
        return '#%d' %  self.anon[cpt_type]

    def _syntax_error(self, string):

        raise ValueError('%s\nExpected format: %s' % (repr(rule)))

    def parse(self, string):
        """Parse string and create object"""

        string = string.strip()
        if string == '':
            return None

        if string[0] in self.comments:
            return None

        self.string = string
        fields = string.split(';')

        fields = self.param_pattern.findall(fields[0])

        # Strip {}, perhaps should do with regexp.
        for m, field in enumerate(fields):
            if field[0] == '{':
                fields[m] = fields[m][1:-1]
        
        name = fields.pop(0)
        match = self.cpt_pattern.match(name)
        if match is None:
            raise ValueError('Unknown component for %s' % name)

        groups = match.groups()
        
        # Add id if anonymous.
        cpt_type, cpt_id = groups[0], groups[1]
        if cpt_id is None:
            name += self._anon_cpt_id(cpt_type)

        # This is the most hackery aspect of this parser where we
        # choose the rule pattern based on a keyword.  If the
        # keyword is not present, default to first rule pattern.
        rule = self.ruledir[cpt_type][0]
        for rule1 in self.ruledir[cpt_type]:
            pos = rule1.pos
            if pos is None:
                continue
            if len(fields) > pos and fields[pos].lower() == rule1.params[pos]:
                rule = rule1
                break

        cpt = rule.process(self.cpts, self.paramdir, name, string, fields)
        
        cpt.id = cpt_id
        cpt.type = cpt_type
        cpt.classname = rule.classname
        return cpt


