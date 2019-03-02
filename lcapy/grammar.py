"""
This module defines a grammar for SPICE-like netlists.

Copyright 2015, 2016 Michael Hayes, UCECE
"""

# SPICE also considers = a delimiter.
delimiters = ' \t(),'

# Comment characters; these must be in the first column.
comments = r'#%*'

# Optional params are in square brackets.
rules = r"""
AM: AMname Np Nm; Ammeter
BAT: BATname Np Nm [Value]; Battery
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
F: Fname Np Nm Vcontrol [Value]; Current controlled current source (note the control current is specified through a voltage source)
FB: FBname Np Nm; Ferrite bead
G: Gname Np Nm Ncp Ncm [Value]; Voltage controlled current source
GY: GYname Np Nm Ncp Ncm [Value]; Gyrator
H: Hname Np Nm Vcontrol [Value]; Current controlled voltage source (note the control current is specified through a voltage source)
I: Iname Np Nm [Value]; Current source
sI: Iname Np Nm s [Value]; s-domain current source
Idc: Iname Np Nm dc [Value]; DC current source
Istep: Iname Np Nm step [Value]; Step current source
Iac: Iname Np Nm ac [Value] [Phase] [Freq]; AC current source
Isin: Iname Np Nm sin Io Ia fo [td] [alpha] [Phase]; Sinusoidal current source
Inoise: Iname Np Nm noise [Value] [NID]; Noise current source
J: Jname Nd Ng Ns [Value]; N channel JFET
Jnjf: Jname Nd Ng Ns njf [Value]; N channel JFET
Jpjf: Jname Nd Ng Ns pjf [Value]; P channel JFET
K: Kname Lname1 Lname2 [Value]; Mutual inductance
L: Lname Np Nm [Value] [IC]; Inductance
M: Mname Nd Ng Ns [Value]; N channel MOSFET
Mnmos: Mname Nd Ng Ns nmos [Value]; N channel MOSFET
Mpmos: Mname Nd Ng Ns pmos [Value]; P channel MOSFET
MT: MTname Np Nm; Motor
MX: MXname P P P; Mixer
O: Oname Np Np [Value]; Open circuit
P: Pname Np Np [Value]; Port
Q: Qname Nc Nb Ne [Value]; NPN transistor
Qnpn: Qname Nc Nb Ne npn [Value]; NPN transistor
Qpnp: Qname Nc Nb Ne pnp [Value]; PNP transistor 
R: Rname Np Nm [Value]; Resistor
Sbox: Sname box; Box
Scircle: Sname circle; Circle
Sellipse: Sname ellipse; Ellipse
Striangle: Sname triangle; Triangle
SPpp: SPname pp P P P; Summing point
SPpm: SPname pm P P P; Summing point
SPppp: SPname ppp P P P P; Summing point
SPpmm: SPname pmm P P P P; Summing point
SPppm: SPname ppm P P P P; Summing point
SWnc: SWname Np Nm nc; Switch normally closed
SWno: SWname Np Nm no; Switch normally open
SWpush: SWname Np Nm push; Pushbutton switch
SWspdt: SWname Nc Np Nm spdt; SPDT switch
TF: TFname Np Nm Ncp Ncm [Value]; Ideal transformer (works to DC!)
TFcore: TFname Np Nm Ncp Ncm core [Value]; Transformer with core (works to DC!)
TFtap: TFname Np Nm Ncp Ncm tap Nt Nt [Value]; Tapped transformer (works to DC!)
TFtapcore: TFname Np Nm Ncp Ncm tapcore Nt Nt [Value]; Tapped transformer with core (works to DC!)
TL: TLname Np Nm Ncp Ncm [Value]; Transmission line
TP: TPname Np Nm Ncp Ncm [Value]; Two port
TR: TRname Pi Po [Value]; Transfer function
Ubuffer: Uname buffer Pi PVss Po PVdd; Buffer with power supplies
Uinverter: Uname inverter Pi PVss Po PVdd; Inverter with power supplies
Udiffamp: Uname diffamp P P P P P; Differential amplifier
Uadc: Uname adc P P P P P P P P; ADC
Udac: Uname dac P P P P P P P P; DAC
Ubox: Uname box P P; Box
Ucircle: Uname circle P P; Circle
Ubox4: Uname box4 P P P P; Box
Ubox12: Uname box12 P P P P P P P P P P P P; Box
Ucircle4: Uname circle4 P P P P; Circle
Uchip1310: Uname chip1310 P P P P P; Chip
Uchip2121: Uname chip2121 P P P P P P; Chip
Uchip3131: Uname chip3131 P P P P P P P P; Chip
Uchip4141: Uname chip4141 P P P P P P P P P P; Chip
V: Vname Np Nm [Value]; Voltage source
sV: Vname Np Nm s [Value]; s-domain voltage source
Vdc: Vname Np Nm dc [Value]; DC voltage source
Vstep: Vname Np Nm step [Value]; Step voltage source
Vac: Vname Np Nm ac [Value] [Phase] [Freq]; AC voltage source
Vsin: Vname Np Nm sin Vo Va fo [td] [alpha] [Phase]; Sinusoidal voltage source
Vnoise: Vname Np Nm noise [Value] [NID]; Noise voltage source
VM: VMname Np Nm; Voltmeter
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
core: keyword;
dc: keyword;
noise: keyword;
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
tap: keyword;
tapcore: keyword;
opamp: keyword;
fdopamp: keyword;
buffer: keyword;
pbuffer: keyword;
pinverter: keyword;
inverter: keyword;
adc: keyword;
dac: keyword;
diffamp: keyword;
box: keyword;
circle: keyword;
ellipse: keyword;
triangle: keyword;
box4: keyword;
box12: keyword;
circle4: keyword;
chip1310: keyword;
chip2121: keyword;
chip3131: keyword;
chip4141: keyword;
pp: keyword;
pm: keyword;
ppp: keyword;
pmm: keyword;
ppm: keyword;
push: keyword;
P: pin; Pin
Pi: pin; Input pin
Po: pin; Output pin
PVdd: pin; Vdd pin
PVss: pin; Vss pin
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
Nt: node; Tap node
Nodelist: nodelist; List of nodes
Phase: value; AC phase
Freq: value; AC frequency
Vo: value; DC voltage offset
Va: value; Sinewave voltage amplitude
Io: value; DC current offset
Ia: value; Sinewave current amplitude
fo: value; Sinewave frequency
td: value; Time delay
alpha: value; Damping factor
Value: value; Value
IC: value; Initial condition
NID: value; Noise identifier
Lname1: name; Inductor1 name 
Lname2: name; Inductor2 name 
Vcontrol: name; Control voltage name 
"""
