import re

sub_super_pattern = re.compile(r"([_\^]){([a-zA-Z]+)([0-9]*)}")


class Latex(object):

    words = ('in', 'out', 'ref', 'rms', 'load', 'source', 'avg',
             'mean', 'peak', 'pp', 'min', 'max', 'src', 'bat',
             'cc', 'ee', 'dd', 'ss', 'ih', 'il', 'oh', 'ol')

    def __init__(self, string):

        self.str = string

    def mathrm(self):
        """Place words in sub- or super-scripts inside a mathrm.
        For example V_{rms} -> V_{\mathrm{rms}}"""
        
        
        def foo(match):
            
            word = match.group(2)
            suffix = word + match.group(3)
            
            if word.lower() in self.words:
                suffix = r'{\mathrm{%s}}' % suffix
            else:
                suffix = r'{%s}' % suffix
            return match.group(1) + suffix

        return sub_super_pattern.sub(foo, self.str)

    def __str__(self):
        
        return self.mathrm()


def latex_str(string):

    return Latex(string).__str__()


def format_label(s):

    if s == '':
        return s

    # If have $ in string then assume necessary parts are in math-mode.
    if '$' in s:
        return s

    # If have _, ^, or \ need to be in math-mode.
    if '_' in s or '^' in s or '\\' in s:
        return '$' + latex_str(s) + '$'
    return s

