"""
This module defines a grammar for SPICE-like netlists.

Copyright 2015, 2016 Michael Hayes, UCECE
"""

# SPICE also considers = a delimiter.
delimiters = r' \t\(\),'

# Optional args are in square brackets.
rules = r"""
D: Dname Np Nm; Diode
Dled: Dname Np Nm led; Light emitting diode
Dzener: Dname Np Nm zener; Zener diode
Dphoto: Dname Np Nm photo; Photo diode
Dtunnel: Dname Np Nm tunnel; Tunnel diode
Dshottky: Dname Np Nm shottky; Shottky diode
E: Ename Np Nm Ncp Ncm [Value]; Voltage controlled voltage source
F: Fname Np Nm Vcontrol [Value]; Voltage controlled current source
G: Gname Np Nm Ncp Ncm [Value]; Current controlled voltage source
H: Hname Np Nm Vcontrol [Value]; Current controlled current source
I: Iname Np Nm [Value]; Current source
Idc: Iname Np Nm dc [Value]; DC current source
Iac: Iname Np Nm ac [Value] [Phase]; AC current source
Isin: Iname Np Nm sin Io Ia fo [td] [alpha] [Phase]; Sinusoidal current source
J: Jname Nd Ng Ns; N channel JFET
Jnjf: Jname Nd Ng Ns njf; N channel JFET
Jpjf: Jname Nd Ng Ns pjf; P channel JFET
K: Kname Lname1 Lname2 [Value]; Mutual inductance
L: Lname Np Nm [Value]; Inductance
M: Mname Nd Ng Ns; N channel MOSFET
Mnmos: Mname Nd Ng Ns nmos; N channel MOSFET
Mpmos: Mname Nd Ng Ns pmos; P channel MOSFET
Q: Qname Nc Nb Ne; NPN transistor
Qnpn: Qname Nc Nb Ne npn; NPN transistor
Qpnp: Qname Nc Nb Ne pnp; PNP transistor 
R: Rname Np Nm [Value]; Resistor
SWnc: SWname Np Nm nc; Switch normally closed
SWno: SWname Np Nm no; Switch normally open
SWpush: SWname Np Nm push; Pushbutton switch
V: Vname Np Nm [Value]; Voltage source
Vdc: Vname Np Nm dc [Value]; DC voltage source
Vac: Vname Np Nm ac [Value] [Phase]; AC voltage source
Vsin: Vname Np Nm sin Vo Va fo [td] [alpha] [Phase]; Sinusoidal voltage source
"""

args = r"""
led: keyword;
zener: keyword;
photo: keyword;
tunnel: keyword;
shottky: keyword;
ac: keyword;
dc: keyword;
sin: keyword;
njf: keyword;
pjf: keyword;
npn: keyword;
pnp: keyword;
nmos: keyword;
pmos: keyword;
no: keyword;
nc: keyword;
push: keyword;
Nb: node; Base node
Nc: node; Collector node
Ncp: node; Positive control node
Ncm: node; Negative control node
Nd: node; Drain node
Ne: node; Emitter node
Ng: node; Gate node
Nm: node; Negative node
Np: node; Positive node
Ns: node; Source node
Phase: value; AC Phase
Vo: value; DC voltage offset
Va: value; Sinewave voltage amplitude
Io: value; DC current offset
Ia: value; Sinewave current amplitude
fo: value; Sinewave frequency
td: value; Time delay
alpha: value; Damping factor
Value: value; Value
Lname1: name; Inductor1 name 
Lname2: name; Inductor2 name 
Vcontrol: name; Control voltage name 
"""
