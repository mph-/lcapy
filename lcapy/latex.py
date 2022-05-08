"""This module provides the Latex class for formatting
values and labels.

Copyright 2020--2022 Michael Hayes, UCECE

"""

import re
from .config import subscripts, greek_letter_names

sub_super_pattern = re.compile(r"([_\^]){([a-zA-Z]+)([0-9]*)}")
greek_letter_name_pattern = re.compile(
    r"([^a-zA-Z\\]|^)(%s)([^a-zA-Z]|$)" % '|'.join(greek_letter_names))
math_symbols_pattern = re.compile(
    '|'.join(['_',  r'\^', r'\\left', r'\\math', r'\\frac', r'\\sqrt']))


class Latex(object):

    def __init__(self, string):

        self.str = string

    def mathrm(self, s):
        """Place words in sub- or super-scripts inside a mathrm.
        For example: V_{rms} -> V_{\mathrm{rms}}"""

        def foo(match):

            word = match.group(2)
            suffix = word + match.group(3)

            # Perhaps look up dictionary to find valid words?
            # Assume that if length 3 or more then a word.
            if word.lower() in subscripts or len(word) > 2:
                suffix = r'{\mathrm{%s}}' % suffix
            else:
                suffix = r'{%s}' % suffix
            return match.group(1) + suffix

        return sub_super_pattern.sub(foo, s)

    def greek(self, s):
        """Replace alpha with \alpha, etc."""

        def foo(match):
            return match.group(1) + '\\' + match.group(2) + match.group(3)

        return greek_letter_name_pattern.sub(foo, s)

    def __str__(self):

        s = self.str
        s = self.greek(s)
        return self.mathrm(s)


def latex_str(string):

    return Latex(string).__str__()


def latex_format_label(s):

    if s == '':
        return s

    # With leading $ and no trailing $, e.g., v=$V1, try to
    # automagically convert to LaTeX string, otherwise pass
    # verbatim.  Currently, this only converts words in sub- or
    # super- scripts to roman.  TODO, handle more cases.
    if s[0] == '$' and s[-1] != '$':
        return '$' + latex_str(s[1:]) + '$'

    # If find a $ assume that the label is correctly formatted.
    if '$' in s:
        return s

    s = latex_str(s)

    s2 = s
    for macro in (r'\hspace', r'\vspace', r'\ ', r'\,', r'\;', r'\!',
                  r'\tiny', r'\scriptsize', r'\footnotesize', r'\small',
                  r'\normalsize', r'\large', '\Large',
                  r'\LARGE', r'\huge', r'\Huge'):
        s2 = s2.replace(macro, '')

    # If have _, ^, \frac, etc.  need to be in math-mode.  Should
    # prevent lcapy generating such strings and warn user to explicity
    # use math mode.  The tricky part is that arbitrary signals may
    # have math symbols, say sqrt.
    math_symbols = ('_',  '^', '\\')
    for math_symbol in math_symbols:
        if math_symbol in s2:
            return '$' + s + '$'
    return s


def latex_format_node_label(s):

    if s[0] == '_':
        return s.replace('_', r'\_')

    parts = s.split('_')
    if len(parts) > 2:
        raise ValueError('Double subscript in %s' % s)
    elif len(parts) == 2:
        return '%s$_{\mathrm{%s}}$' % (parts[0], parts[1])

    return s
