from nose.tools import *
import sys
sys.path.append('..')

from lcapy.parser import Parser
import lcapy.schemcpts as schemcpts
import lcapy.grammar as grammar

parser = Parser(schemcpts, grammar)
parse = parser.parse

@raises(ValueError)
def test_Exception1():
    '''Test missing arg'''

    parse('V1 2')


@raises(ValueError)
def test_Exception2():
    '''Test unknown component'''

    parse('A1 1 2')

    
