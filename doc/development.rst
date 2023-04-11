.. _development:

===========
Development
===========


Debugging expressions
=====================

Lcapy expression have a `pdb()` method that invokes the Python
debugger.  For example, to debug the inverse Laplace transform for the
expression `1 / s` use::

   >>> (1 / s).pdb().ILT()


Import timing
=============

To find which modules are slow to import, use::

   $ python3 -X importtime -c "import lcapy" 2> /tmp/lcapy.log
   $ tuna /tmp/lcapy.log

This requires Python3.7.  tuna graphically shows which modules take the most time to import.   tuna can be installed with pip::

  $ pip3 install tuna


Adding new components
=====================

There are three steps:
1. Define component in `grammar.py`
2. Add associated class in `mnacpts.py` that handles the electrical properties
3. Add associated class in `schemcpts.py` that draws the component

Note, both `mnacpts.py` and `schemcpts.py` are in dire need of splitting into multiple files; ideally one per class.


Python package dependencies
===========================

To find Python package dependencies::

   $ pip3 install pipdeptree
   $ pipdeptree
