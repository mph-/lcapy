from mcircuit import Opamp, R, Series

# Create simple inverting amplifier
Rf = 1000
Ri = 50

a = Opamp()

# Connect V+ to ground.
b = a.shortcircuit(1)

# Add feedback resistor.
c = b.bridge(R(Rf))

# Add input resistor to V-.
d = c.prepend(Series(R(Ri)))

print d.Vtransfer

