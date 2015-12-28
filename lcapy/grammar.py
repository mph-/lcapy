grammar = r"""
start: c | d | dled | dphoto | dzener | dshottky | dtunnel | e | f | g | h | i | idc | iac | isin | jnjf | jpjf | k | l | qnpn | qpnp | r | swnc | swno | swpush | v | vdc | vac | vsin;
c: cname pnode nnode value? opts?;
cname: 'C(\d|\w)+';
d: dname pnode nnode opts?;
dled: dname pnode nnode 'led' opts?;
dzener: dname pnode nnode 'zener' opts?;
dphoto: dname pnode nnode 'photo' opts?;
dtunnel: dname pnode nnode 'tunnel' opts?;
dshottky: dname pnode nnode 'shottky' opts?;
dname: 'D(\d|\w)+';
e: ename pnode nnode cpnode cnnode value? opts?;
ename: 'E(\d|\w)+';
f: fname pnode nnode vname value? opts?;
fname: 'F(\d|\w)+';
g: gname pnode nnode cpnode cnnode value? opts?;
gname: 'G(\d|\w)+';
h: hname pnode nnode vname value? opts?;
hname: 'H(\d|\w)+';
i: iname pnode nnode value? opts?;
idc: iname pnode nnode 'dc' value? opts?;
iac: iname pnode nnode 'ac' value? phase? opts?;
isin: iname pnode nnode 'sin' io ia fo td? alpha? phase? opts?;
io: val;
ia: val;
iname: 'I(\d|\w)+';
jnjf: jname dnode gnode snode 'njf'? opts?;
jpjf: jname dnode gnode snode 'pjf' opts?;
jname: 'J(\d|\w)+';
k: kname lname lname value? opts?;
kname: 'K(\d|\w)+';
l: lname pnode nnode value? opts?;
lname: 'L(\d|\w)+';
// C B E or D G S
qnpn: qname cnode bnode enode 'npn'? opts?;
qpnp: qname cnode bnode enode 'pnp' opts?;
qname: 'Q(\d|\w)+';
r: rname pnode nnode value? opts?;
rname: 'R(\d|\w)+';
swnc: swname pnode nnode 'nc' opts?;
swno: swname pnode nnode 'no' opts?;
swpush: swname pnode nnode 'push' opts?;
swname: 'SW(\d|\w)+';
v: vname pnode nnode value? opts?;
vdc: vname pnode nnode 'dc' value? opts?;
vac: vname pnode nnode 'ac' value? phase? opts?;
vsin: vname pnode nnode 'sin' vo va fo td? alpha? phase? opts?;
vo: val;
va: val;
fo: val;
td: val;
alpha: val;
phase: val;
// opts: ;
keyname: 'dir' | 'colour' | 'size' | 'v' | 'l' | 'i';
keyword: 'mirror' | 'reverse' | 'down' | 'up' | 'left' | 'right' | 'implict' | 'ground' | 'sground' ;
opts: ';' opt opt* ;
keyval: '=.*';
keypair: keyname keyval ;
@opt: keyword | keypair;
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
//WHITESPACE: '[ \t\(\)=,]+' (%ignore);
WHITESPACE: '[ \t\(\),]+' (%ignore);
"""
