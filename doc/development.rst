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


Python package dependencies
===========================

To find Python package dependencies::

   $ pip3 install pipdeptree
   $ pipdeptree
