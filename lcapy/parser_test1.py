from parser import Parser
import schemcpts as cpts

parser = Parser(cpts)
parse = parser.parse

r = parse('R1 2 3 4')
print(r)

r1 = parse('R1 2 3')
print(r1)

q = parse('Q1 2 3 4 pnp')
print(q)

j = parse('Ja2 2 3 4 pjf')
print(j)

s = parse('SW1 2 3 push')
print(s)

#d1 = parse('D1 2 3 led; dir : down')
d1 = parse('D1 2 3 led; dir=down')
print(d1)

#v1 = parse('V1 2 3 dc=4')
v1 = parse('V1 2 3 dc 4')
print(v1)

v2 = parse('V1 2 3 ac 10, 20')
print(v2)

v3 = parse('V1 2 3 sin(0 10 1000)')
print(v3)

v4 = parse('V1 2 3 1.0e3')
print(v4)

v5 = parse('V1 2 3 1e3')
print(v5)

v6 = parse('V1 2 3; right colour=blue size=1.5 l=$V_1$')
print(v6)

r2 = parse('R1 2 3 {5 * a}')
print(r2)

i3 = parse('I1 2 3 sin(0 10 1000)')
print(i3)

e1 = parse('E1 2 3 4 5 1e3')
print(e1)
