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
    '''Test too many args'''

    parse('V1 2 3 4 5 6')

@raises(ValueError)
def test_Exception3():
    '''Test too many args'''

    parse('V1 2 3 dc 4 5 6')

@raises(ValueError)
def test_Exception4():
    '''Test unknown component'''

    parse('A1 1 2')

def test_V1():
    '''Test voltage source'''
    
    assert_equals(type(parse('V1 1 2')), schemcpts.classes['V'], 'Class not V')

def test_Vdc1():
    '''Test dc voltage source'''
    
    assert_equals(type(parse('V1 1 2 dc')), schemcpts.classes['Vdc'], 'Class not Vdc')

def test_Vac1():
    '''Test ac voltage source'''
    
    assert_equals(type(parse('V1 1 2 ac')), schemcpts.classes['Vac'], 'Class not Vac')

def test_Vsin1():
    '''Test sin voltage source'''
    
    assert_equals(type(parse('V1 1 2 sin(1, 2, 3)')), schemcpts.classes['Vsin'], 'Class not Vsin')

def test_Vexpr():
    '''Test voltage source with expr'''
    
    assert_equals(type(parse('V1 1 2 {a * 5}')), schemcpts.classes['V'], 'Class not V')

def test_Vac2():
    '''Test ac voltage source'''
    
    assert_equals(type(parse('V1 1 2 ac {1j} {0} 3')), schemcpts.classes['Vac'], 'Class not Vac')    

def test_Vquotes():
    '''Test voltage source with arg in quotes'''
    
    assert_equals(type(parse('V1 1 2 "a * 5"')), schemcpts.classes['V'], 'Class not V')

# def test_opamp():
#     '''Test opamp'''
#    
#     assert_equals(type(parse('E 1 2 opamp 3 4')), schemcpts.classes['Eopamp'], 'Class not Eopamp')
