grammar = r"""
start: c | d | dled | dphoto | dzener | dshottky | dtunnel | e | f | g | h | i | idc | iac | isin | jnjf | jpjf | k | l | qnpn | qpnp | r | swnc | swno | swpush | v | vdc | vac | vsin;
c: cname pnode nnode value?;
cname: 'C(\d|\w)+';
d: dname pnode nnode;
dled: dname pnode nnode 'led';
dzener: dname pnode nnode 'zener';
dphoto: dname pnode nnode 'photo';
dtunnel: dname pnode nnode 'tunnel';
dshottky: dname pnode nnode 'shottky';
dname: 'D(\d|\w)+';
e: ename pnode nnode cpnode cnnode value?;
ename: 'E(\d|\w)+';
f: fname pnode nnode vname value?;
fname: 'F(\d|\w)+';
g: gname pnode nnode cpnode cnnode value?;
gname: 'G(\d|\w)+';
h: hname pnode nnode vname value?;
hname: 'H(\d|\w)+';
i: iname pnode nnode value?;
idc: iname pnode nnode 'dc' value?;
iac: iname pnode nnode 'ac' value? phase?;
isin: iname pnode nnode 'sin' io ia fo td? alpha? phase?;
io: val;
ia: val;
iname: 'I(\d|\w)+';
jnjf: jname dnode gnode snode 'njf'?;
jpjf: jname dnode gnode snode 'pjf';
jname: 'J(\d|\w)+';
k: kname lname lname value?;
kname: 'K(\d|\w)+';
l: lname pnode nnode value?;
lname: 'L(\d|\w)+';
// C B E or D G S
qnpn: qname cnode bnode enode 'npn'?;
qpnp: qname cnode bnode enode 'pnp';
qname: 'Q(\d|\w)+';
r: rname pnode nnode value?;
rname: 'R(\d|\w)+';
swnc: swname pnode nnode 'nc';
swno: swname pnode nnode 'no';
swpush: swname pnode nnode 'push';
swname: 'SW(\d|\w)+';
v: vname pnode nnode value?;
vdc: vname pnode nnode 'dc' value?;
vac: vname pnode nnode 'ac' value? phase?;
vsin: vname pnode nnode 'sin' vo va fo td? alpha? phase?;
vo: val;
va: val;
fo: val;
td: val;
alpha: val;
phase: val;
vname: 'V(\d|\w)+';
@integer: '\d+';
@float: '-?([1-9]\d*|\d)\.(\d+)?([eE][+-]?\d+)?' | '-?([1-9]\d*|\d)[eE]([+-]?\d+)?';
value: float | integer | '\{.*\}';
@val: float | integer | '\{.*\}';
// Positive node
pnode: xnode;
// Negative node
nnode: xnode;
// Positive controlling node
cpnode: xnode;
// Negative controlling node
cnnode: xnode;
// Collector node
cnode: xnode;
// Base node
bnode: xnode;
// Emitter node
enode: xnode;
// Drain node
dnode: xnode;
// Gate node
gnode: xnode;
// Source node
snode: xnode;
@xnode: '\d+';
@nul: ;
WHITESPACE: '[ \t\(\)=,]+' (%ignore);
"""
