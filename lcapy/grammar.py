"""
This module defines a grammar for SPICE-like netlists.

Copyright 2015--2022 Michael Hayes, UCECE
"""

# SPICE also considers = a delimiter.
delimiters = ' \t(),'

# Comment characters; these must be in the first column.
comments = r'#%*'

# Each line defines a rule.  The field before the colon is the classname.  This
# must be unique.
# The second field (up to name) is used to match the rule along with the optional
# keyword parameter.
# Optional parameters are in square brackets.
rules = r"""
A: Aname Np; Annotation
ADC: ADCname Np Nm; ADC
AM: AMname Np Nm; Ammeter
BAT: BATname Np Nm [Value=name]; Battery
BL: BLname Np Nm; Block
C: Cname Np Nm [Value=name] [IC]; Capacitor
Cable: Cablename; Cable
CPE: CPEname Np Nm [Value=name] [Power=1]; Constant phase element
D: Dname Np Nm; Diode
DAC: DACname Np Nm; DAC
Dled: Dname Np Nm led; Light emitting diode
Dzener: Dname Np Nm zener; Zener diode
Dphoto: Dname Np Nm photo; Photo diode
Dtunnel: Dname Np Nm tunnel; Tunnel diode
Dschottky: Dname Np Nm schottky; Schottky diode
E: Ename Np Nm Ncp Ncm [Value=name] [Ac=0]; Voltage controlled voltage source
VCVS: VCVSname Np Nm Ncp Ncm [Value=name]; Voltage controlled voltage source
Eopamp: Ename Np Nm opamp Ncp Ncm [Ad=name] [Ac=0] [Ro=0]; Opamp
Efdopamp: Ename Np Nm fdopamp Ncp Ncm Nocm [Ad=name] [Ac=0]; Fully differential opamp
Einamp: Ename Np Nm inamp Ncp Ncm NRp NRm [Ad=name] [Ac=0] [Rf=Rf]; Instrumentation opamp
Eamp: Ename Np Nm amp Ncp Ncm [Ad=name] [Ac=0]; Amplifier
F: Fname Np Nm Vcontrol [Value=name]; Current controlled current source (note the control current is specified through a voltage source)
CCCS: CCCSname Np Nm Vcontrol [Value=name]; Current controlled current source (note the control current is specified through a voltage source)
FB: FBname Np Nm; Ferrite bead
FS: FSname Np Nm; Fuse
G: Gname Np Nm Ncp Ncm [Value=name]; Voltage controlled current source
VCCS: VCCSname Np Nm Ncp Ncm [Value=name]; Voltage controlled current source
GY: GYname Np Nm Ncp Ncm [Value=name]; Gyrator
H: Hname Np Nm Vcontrol [Value=name]; Current controlled voltage source (note the control current is specified through a voltage source)
CCVS: CCVSname Np Nm Vcontrol [Value=name]; Current controlled voltage source (note the control current is specified through a voltage source)
I: Iname Np Nm [Value=name]; Current source
sI: Iname Np Nm s [Value=name]; s-domain current source
Idc: Iname Np Nm dc [Value=name]; DC current source
Istep: Iname Np Nm step [Value=name]; Step current source
Iac: Iname Np Nm ac [Value=name] [Phase] [Freq]; AC current source
Isin: Iname Np Nm sin Io Ia fo [td] [alpha] [Phase]; Sinusoidal current source
Inoise: Iname Np Nm noise [Value=name] [NID]; Noise current source
J: Jname Nd Ng Ns [Value=name]; N channel JFET
Jnjf: Jname Nd Ng Ns njf [Value=name]; N channel JFET
Jpjf: Jname Nd Ng Ns pjf [Value=name]; P channel JFET
k: kname Np Nm [Value=name] [IC]; Spring
K: Kname Lname1 Lname2 [Value=name]; Mutual inductance (L1 on left, L2 on right)
L: Lname Np Nm [Value=name] [IC]; Inductance
m: mname Np Nm [Value=name] [IC]; Mass
M: Mname Nd Ng Ns [Value=name]; N channel MOSFET
Mnmos: Mname Nd Ng Ns nmos [Value=name]; N channel MOSFET
Mpmos: Mname Nd Ng Ns pmos [Value=name]; P channel MOSFET
MISC: MISCname Np Nm; Miscellaneous circuitikz bipole
MT: MTname Np Nm; Motor
MX: MXname P P P; Mixer
NR: NRname Np Nm [Value=name]; Noiseless resistor
O: Oname Np Np; Open circuit
P: Pname Np Np; Port
Q: Qname Nc Nb Ne [Value=name]; NPN transistor
Qnpn: Qname Nc Nb Ne npn [Value=name]; NPN transistor
Qpnp: Qname Nc Nb Ne pnp [Value=name]; PNP transistor
r: rname Np Nm [Value=name]; Damper
R: Rname Np Nm [Value=name]; Resistor
REL: RELname Np Nm [Value=name]; Reluctance
RV: RVname Np Nm No [Value=name] [Value=name]; Potentiometer
Sbox: Sname box; Box
Scircle: Sname circle; Circle
Sellipse: Sname ellipse; Ellipse
Striangle: Sname triangle; Triangle
SPpp: SPname pp P P P; Summing point
SPpm: SPname pm P P P; Summing point
SPppp: SPname ppp P P P P; Summing point
SPpmm: SPname pmm P P P P; Summing point
SPppm: SPname ppm P P P P; Summing point
SW: SWname Np Nm [Time=0]; Switch normally open
SWnc: SWname Np Nm nc [Time=0]; Switch normally closed
SWno: SWname Np Nm no [Time=0]; Switch normally open
SWpush: SWname Np Nm push [Time=0]; Pushbutton switch
SWspdt: SWname Nc Np Nm spdt [Time=0]; SPDT switch
TF: TFname Np Nm Ncp Ncm [Value=name]; Ideal transformer (works to DC!)
TFcore: TFname Np Nm Ncp Ncm core [Value=name]; Transformer with core (works to DC!)
TFtap: TFname Np Nm Ncp Ncm tap Nt Nt [Value=name]; Tapped transformer (works to DC!)
TFtapcore: TFname Np Nm Ncp Ncm tapcore Nt Nt [Value=name]; Tapped transformer with core (works to DC!)
TL: TLname Np Nm Ncp Ncm [Z0=Z0] [Gamma=Gamma] [Length=l]; Transmission line
TLlossless: TLname Np Nm Ncp Ncm lossless [Z0=Z0] [Gamma=Gamma] [Length=l]; Lossless transmission line
TP: TPname Np Nm Ncp Ncm; Generic two-port
TPA: TPname Np Nm Ncp Ncm A A11 A12 A21 A22 [V1] [I1]; A-parameter two-port
TPB: TPname Np Nm Ncp Ncm B B11 B12 B21 B22 [V2] [I2]; B-parameter two-port
TPG: TPname Np Nm Ncp Ncm G G11 G12 G21 G22 [I1] [V2]; G-parameter two-port
TPH: TPname Np Nm Ncp Ncm H H11 H12 H21 H22 [V1] [I2]; H-parameter two-port
TPY: TPname Np Nm Ncp Ncm Y Y11 Y12 Y21 Y22 [I1] [I2]; Y-parameter two-port
TPZ: TPname Np Nm Ncp Ncm Z Z11 Z12 Z21 Z22 [V1] [V2]; Z-parameter two-port
TR: TRname Pi Po [Value=name]; Transfer function
TVtriode: TVname Nanode Ngrid Ncathode [mu] [gm] [Cgk] [Cga]; Thermionic valve triode
Ubuffer: Uname buffer; Buffer
Uinverter: Uname inverter; Inverter
Udiffamp: Uname diffamp; Differential amplifier
Udiffdriver: Uname diffdriver; Differential driver
Uopamp: Uname opamp; Opamp
Uinamp: Uname inamp; Instrumentation amplifier
Uisoamp: Uname isoamp; Isolated amplifier
Ufdopamp: Uname fdopamp; Fully differential opamp
Uregulator: Uname regulator; Voltage regulator
Uadc: Uname adc; ADC
Udac: Uname dac; DAC
Ubox: Uname box; Box
Ucircle: Uname circle; Circle
Ubox4: Uname box4; Box
Ubox12: Uname box12; Box
Ucircle4: Uname circle4; Circle
Uchip1313: Uname chip1313; Chip
Uchip2121: Uname chip2121; Chip
Uchip2222: Uname chip2222; Chip
Uchip3131: Uname chip3131; Chip
Uchip3333: Uname chip3333; Chip
Uchip4141: Uname chip4141; Chip
Uchip4444: Uname chip4444; Chip
Uchip8181: Uname chip8181; Chip
Uchip8888: Uname chip8888; Chip
Umux21: Uname mux21; 2-1 multiplexer
Umux41: Uname mux41; 4-1 multiplexer
Umux42: Uname mux42; 4-2 multiplexer
Udff: Uname dff; D flip/flop
Ujkff: Uname jkff; D flip/flop
Urslatch: Uname rslatch; RS latch
V: Vname Np Nm [Value=name]; Voltage source
sV: Vname Np Nm s [Value=name]; s-domain voltage source
Vdc: Vname Np Nm dc [Value=name]; DC voltage source
Vstep: Vname Np Nm step [Value=name]; Step voltage source
Vac: Vname Np Nm ac [Value=name] [Phase] [Freq]; AC voltage source
Vsin: Vname Np Nm sin Vo Va fo [td] [alpha] [Phase]; Sinusoidal voltage source
Vnoise: Vname Np Nm noise [Value=name] [NID]; Noise voltage source
VM: VMname Np Nm; Voltmeter
W: Wname Np Np; Wire
XT: XTname Np Nm; Crystal
Y: Yname Np Np [Value=name]; Admittance
Z: Zname Np Np [Value=name]; Impedance
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
lossless: keyword;
opamp: keyword;
inamp: keyword;
isoamp: keyword;
fdopamp: keyword;
amp: keyword;
regulator: keyword;
buffer: keyword;
pbuffer: keyword;
pinverter: keyword;
inverter: keyword;
diffdriver: keyword;
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
chip1313: keyword;
chip2121: keyword;
chip2222: keyword;
chip3131: keyword;
chip3333: keyword;
chip4141: keyword;
chip4444: keyword;
chip8181: keyword;
chip8888: keyword;
mux21: keyword;
mux41: keyword;
mux42: keyword;
pp: keyword;
pm: keyword;
ppp: keyword;
pmm: keyword;
ppm: keyword;
push: keyword;
dff: keyword; D flip-flop
jkff: keyword; JK flip-flop
rslatch: keyword; RS latch
A: keyword;
B: keyword;
G: keyword;
H: keyword;
Y: keyword;
Z: keyword;
P: pin; Pin
Pi: pin; Input pin
Po: pin; Output pin
Nanode: node; Anode node
Ngrid: node; Grid node
Ncathode: node; Cathode node
mu: value; Amplification factor
gm: value; Transconductance
Cgk: value; Grid to cathode capacitance
Cga: value; Grid to anode capacitance
Nb: node; Base node
Nc: node; Collector node
Ncp: node; Positive control node
Ncm: node; Negative control node
Nocm: node; Output common-mode node
Nd: node; Drain node
Ne: node; Emitter node
Ng: node; Gate node
Ni: node; Input node
Nm: node; Negative node
No: node; Output node
Np: node; Positive node
NRp: node; Gain resistor positive node
NRm: node; Gain resistor negative node
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
Time: value; Time
Value: value; Value
IC: value; Initial condition
NID: value; Noise identifier
Power: value; Power
Lname1: name; Inductor1 name
Lname2: name; Inductor2 name
Vcontrol: name; Control voltage name
Z0: value; Characteristic impedance
Gamma: value; Propagation constant
Length: value; Transmission line length
Ac: value; Common-mode gain
Ad: value; Differential gain
Rf: value; Feedback resistance
Ro: value; Output resistance
A11: value; A11
A12: value; A12
A21: value; A21
A22: value; A22
B11: value; B11
B12: value; B12
B21: value; B21
B22: value; B22
G11: value; G11
G12: value; G12
G21: value; G21
G22: value; G22
H11: value; H11
H12: value; H12
H21: value; H21
H22: value; H22
Y11: value; Y11
Y12: value; Y12
Y21: value; Y21
Y22: value; Y22
Z11: value; Z11
Z12: value; Z12
Z21: value; Z21
Z22: value; Z22
V1: value; V1
I1: value; I1
V2: value; V2
I2: value; I2
"""
