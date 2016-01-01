import sys
sys.path.append('../..')

from parser import Parser
import schemcpts
import grammar

parser = Parser(schemcpts, grammar)
parse = parser.parse

r = parse('R1 2 3 4')
print(r)

r1 = parse('R1 2 3')
print(r1)

q1 = parse('Q1 2 3 4 pnp')
print(q1)

j = parse('Ja2 2 3 4 pjf')
print(j)

s = parse('SW1 2 3 push')
print(s)

d1 = parse('D1 2 3 led')
print(d1)

v1 = parse('V1 2 3')
print(v1)

#v1 = parse('V1 2 3 dc=4')
#print(v1)

v1 = parse('V1 2 3 dc 4')
print(v1)

v2 = parse('V2 2 3 ac 10, 20')
print(v2)

v3 = parse('V3 2 3 sin(0 10 1000)')
print(v3)

v4 = parse('V4 2 3 1.0e3')
print(v4)

v5 = parse('V5 2 3 1e3')
print(v5)

v6 = parse('V 2 3 1e3')
print(v6)

r2 = parse('R2 2 3 {5 * a}')
print(r2)

i3 = parse('I3 2 3 sin(0 10 1000)')
print(i3)

r2 = parse('R2 2 3 1e3; dir=down')
print(r2)
