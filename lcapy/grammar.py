"""
This module defines a grammar for SPICE-like netlists.

Copyright 2015, 2016 Michael Hayes, UCECE
"""

# SPICE also considers = a delimiter.
delimiters = r' \t\(\),'

# Comment characters
comments = r'#%'

# Optional params are in square brackets.
rules = r"""
C: Cname Np Nm [Value] [IC]; Capacitor
D: Dname Np Nm; Diode
Dled: Dname Np Nm led; Light emitting diode
Dzener: Dname Np Nm zener; Zener diode
Dphoto: Dname Np Nm photo; Photo diode
Dtunnel: Dname Np Nm tunnel; Tunnel diode
Dschottky: Dname Np Nm schottky; Schottky diode
E: Ename Np Nm Ncp Ncm [Value]; Voltage controlled voltage source
Eopamp: Ename Np Nm opamp Ncp Ncm [Value]; Opamp
Efdopamp: Ename Np Nm fdopamp Ncp Ncm [Value]; Fully differential opamp
F: Fname Np Nm Vcontrol [Value]; Current controlled current source
G: Gname Np Nm Ncp Ncm [Value]; Voltage controlled current source
H: Hname Np Nm Vcontrol [Value]; Current controlled voltage source
I: Iname Np Nm [Value]; Current source
sI: Iname Np Nm s [Value]; s-domain current source
Idc: Iname Np Nm dc [Value]; DC current source
Iac: Iname Np Nm ac [Value] [Phase]; AC current source
Isin: Iname Np Nm sin Io Ia fo [td] [alpha] [Phase]; Sinusoidal current source
J: Jname Nd Ng Ns [Value]; N channel JFET
Jnjf: Jname Nd Ng Ns njf [Value]; N channel JFET
Jpjf: Jname Nd Ng Ns pjf [Value]; P channel JFET
K: Kname Lname1 Lname2 [Value]; Mutual inductance
L: Lname Np Nm [Value] [IC]; Inductance
M: Mname Nd Ng Ns [Value]; N channel MOSFET
Mnmos: Mname Nd Ng Ns nmos [Value]; N channel MOSFET
Mpmos: Mname Nd Ng Ns pmos [Value]; P channel MOSFET
O: Oname Np Np [Value]; Open circuit
P: Pname Np Np [Value]; Port
Q: Qname Nc Nb Ne [Value]; NPN transistor
Qnpn: Qname Nc Nb Ne npn [Value]; NPN transistor
Qpnp: Qname Nc Nb Ne pnp [Value]; PNP transistor 
R: Rname Np Nm [Value]; Resistor
SWnc: SWname Np Nm nc; Switch normally closed
SWno: SWname Np Nm no; Switch normally open
SWpush: SWname Np Nm push; Pushbutton switch
SWspdt: SWname Nc Np Nm spdt; SPDT switch
TF: TFname Np Nm Ncp Ncm [Value]; Transformer (works to DC!)
TL: TLname Np Nm Ncp Ncm [Value]; Transmission line
TP: TPname Np Nm Ncp Ncm [Value]; Two port
V: Vname Np Nm [Value]; Voltage source
sV: Vname Np Nm s [Value]; s-domain voltage source
Ubuffer: Uname Ni No buffer; Buffer
Upbuffer: Uname Ni No Np Nm pbuffer; Buffer with power supplies
Uinverter: Uname Ni No inverter; Inverter
Vdc: Vname Np Nm dc [Value]; DC voltage source
Vstep: Vname Np Nm step [Value]; Step voltage source
Vac: Vname Np Nm ac [Value] [Phase]; AC voltage source
Vsin: Vname Np Nm sin Vo Va fo [td] [alpha] [Phase]; Sinusoidal voltage source
W: Wname Np Np; Wire
XT: XTname Np Nm; Crystal
Y: Yname Np Np [Value]; Admittance
Z: Zname Np Np [Value]; Impedance
"""

params = r"""
led: keyword;
zener: keyword;
photo: keyword;
tunnel: keyword;
schottky: keyword;
s: keyword;
ac: keyword;
dc: keyword;
step: keyword;
sin: keyword;
njf: keyword;
pjf: keyword;
npn: keyword;
pnp: keyword;
nmos: keyword;
pmos: keyword;
no: keyword;
nc: keyword;
spdt: keyword;
opamp: keyword;
fdopamp: keyword;
buffer: keyword;
pbuffer: keyword;
inverter: keyword;
push: keyword;
Nb: node; Base node
Nc: node; Collector node
Ncp: node; Positive control node
Ncm: node; Negative control node
Nd: node; Drain node
Ne: node; Emitter node
Ng: node; Gate node
Ni: node; Input node
Nm: node; Negative node
No: node; Output node
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
IC: value; Initial condition
Lname1: name; Inductor1 name 
Lname2: name; Inductor2 name 
Vcontrol: name; Control voltage name 
"""
