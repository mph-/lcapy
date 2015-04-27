===============
Simple circuits
===============



V-R-C circuit (1)
=================

This example plots the transient voltage across a capacitor in a series R-L circuit:

.. image:: examples/netlists/circuit-VRC1.png
   :width: 7cm

.. literalinclude:: examples/netlists/circuit-VRC1-vc.py

.. image:: examples/netlists/circuit-VRC1-vc.png
   :width: 15cm


V-R-C circuit (2)
=================

This example is the same as the previous example but it uses an
alternative method of plotting.

.. literalinclude:: examples/netlists/circuit-VRC2-vc.py

.. image:: examples/netlists/circuit-VRC2-vc.png
   :width: 15cm


V-R-L-C circuit (1)
===================

This example plots the transient voltage across a resistor in a series R-L-C circuit:

.. image:: examples/netlists/circuit-VRLC1.png
   :width: 7cm

.. literalinclude:: examples/netlists/circuit-VRLC1-vr.py

.. image:: examples/netlists/circuit-VRLC1-vr.png
   :width: 15cm



V-R-L-C circuit (2)
===================

This is the same as the previous example but with a different resistor value giving an underdamped response:

.. image:: examples/netlists/circuit-VRLC2.png
   :width: 7cm

.. literalinclude:: examples/netlists/circuit-VRLC2-vr.py

.. image:: examples/netlists/circuit-VRLC2-vr.png
   :width: 15cm


