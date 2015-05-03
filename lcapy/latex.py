import re

pattern = re.compile(r"([_\^]){([\w]+)}")


class Latex(object):

    def __init__(self, string):

        self.str = string

    def mathrm(self):
        """Place words in sub- or super-scripts inside a mathrm"""
        
        def foo(match):
            
            fred = match.group(2)
            
            if fred in ('in', 'out', 'ref', 'rms', 'load', 'source', 'avg',
                        'mean', 'peak', 'pp', 'min', 'max', 'src'):
                fred = r'{\mathrm{%s}}' % fred
            else:
                fred = r'{%s}' % fred
            return match.group(1) + fred

        return pattern.sub(foo, self.str)

    def __str__(self):
        
        return self.mathrm()


def latex_str(string):

    return Latex(string).__str__()
