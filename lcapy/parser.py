"""This module performs parsing of SPICE-like netlists.  It uses a
custom parser rather than lex/yacc to give better error messages.

Copyright 2015, 2016 Michael Hayes, UCECE

"""

import re

# Could use a script to generate parser and parsing tables if speed
# was important.

class Arg(object):

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

    def __init__(self, cpt_type, classname, args, comment, pos):
        
        self.type = cpt_type
        self.classname = classname
        self.args = args
        self.comment = comment
        self.pos = pos

    def __repr__(self):

        return self.type + 'name ' + ' '.join(self.args)

    def syntax_error(self, error, string):

        raise ValueError('Syntax error: %s when parsing %s\nExpected format: %s' % (error, string, repr(self)))        

    def process(self, cpts, argdir, name, string, fields, opts_string):

        args = self.args
        if len(fields) > len(args):
            self.syntax_error('Too many args', string)

        # Create instance of component object
        try:
            newclass = getattr(cpts, self.classname)
        except:
            newclass = cpts.classes[self.classname]

        cpt = newclass(name, *fields)
        cpt.string = string
        cpt.opts_string = opts_string
        cpt.nodes = ()

        for m, arg in enumerate(args):

            if m >= len(fields):
                # Optional argument
                if arg[0] == '[':
                    break
                self.syntax_error('Missing arg %s' % arg, string)

            if arg[0] == '[':
                arg = arg[1:-1]

            field = fields[m]
            if argdir[arg].base == 'node':
                if field.find('.') != -1:
                    self.syntax_error('Found . in node name %s' % field, string)                
                cpt.nodes += (field, )

            # Perhaps gobble keywords?

            # Add attribute.  Perhaps the __init__ method for
            # the class should create these from the args?
            setattr(cpt, arg, field)

        return cpt


class Parser(object):

    def __init__(self, cpts, grammar):
        """cpts is a module containing a class for each component
        grammar is a module defining the syntax of a netlist"""

        # A string defining the syntax for a netlist
        rules = grammar.rules
        # A string defining arguments
        args = grammar.args
        # A string defining delimiter characters
        delimiters = grammar.delimiters
        # A string defining comment characters
        self.comments = grammar.comments

        self.cpts = cpts
        self._anon_count = 0
        self.argdir = {}
        self.ruledir = {}
        
        for arg in args.split('\n'):
            self._add_arg(arg)

        for rule in rules.split('\n'):
            self._add_rule(rule)

        cpts = self.ruledir.keys()
        cpts.sort(key=len, reverse=True)

        self.cpt_pattern = re.compile('(%s)(\w+)?' % '|'.join(cpts))
        # strings in curly braces are expressions so do not split.
        self.arg_pattern = re.compile('\{.*\}|".*"|[^%s]+' % delimiters)

    def _add_arg(self, string):

        if string == '':
            return

        fields = string.split(':')
        argname = fields[0]
        fields = fields[1].split(';')
        argbase = fields[0].strip()
        comment = fields[1].strip()
        
        self.argdir[argname] = Arg(argname, argbase, comment)

    def _add_rule(self, string):

        if string == '':
            return

        fields = string.split(':')
        cpt_classname = fields[0]
        fields = fields[1].split(';')
        string = fields[0].strip()
        comment = fields[1].strip()
        
        fields = string.split(' ')
        args = fields[1:]

        cpt_type = fields[0][0:-4]
        
        pos = None
        newargs = ()
        for m, arg in enumerate(args):
            if arg[0] == '[':
                arg = arg[1:-1]
            if arg not in self.argdir:
                raise ValueError('Unknown argument %s for %s' % (arg, string))
            if pos is None and self.argdir[arg].base == 'keyword':
                pos = m

        if cpt_type not in self.ruledir:
            self.ruledir[cpt_type] = ()
        self.ruledir[cpt_type] += (Rule(cpt_type, cpt_classname,
                                        args, comment, pos), )

    def _anon_cpt_id(self):

        self._anon_count += 1
        return '#%d' % self._anon_count

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
        opts_string = fields[1].strip() if len(fields) > 1 else ''

        fields = self.arg_pattern.findall(fields[0])
        
        name = fields.pop(0)
        match = self.cpt_pattern.match(name)
        if match is None:
            raise ValueError('Unknown component for %s' % name)

        groups = match.groups()
        
        # Add id if anonymous.
        cpt_type, cpt_id = groups[0], groups[1]
        if cpt_id is None:
            name += self._anon_cpt_id()

        # This is the most hackery aspect of this parser where we
        # choose the rule pattern based on a keyword.  If the
        # keyword is not present, default to first rule pattern.
        rule = self.ruledir[cpt_type][0]
        for rule1 in self.ruledir[cpt_type]:
            pos = rule1.pos
            if pos is None:
                continue
            if len(fields) > pos and fields[pos].lower() == rule1.args[pos]:
                rule = rule1
                break

        cpt = rule.process(self.cpts, self.argdir, name, string,
                           fields, opts_string)
        
        cpt.id = cpt_id
        cpt.type = cpt_type
        return cpt


