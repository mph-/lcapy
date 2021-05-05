
.. lcapy documentation master file, created by
   sphinx-quickstart on Mon Feb 17 17:08:01 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Lcapy's documentation!
=================================

Lcapy is a Python package for linear circuit analysis.  It uses
SymPy (symbolic Python) for symbolic analysis.  As well as circuit analysis, Lcapy can semi-automate the drawing of high-quality schematics from a netlist, including diodes, transistors, and other non-linear components.  It can also symbolically transform both continuous-time and discrete-time signals, using Laplace, Fourier, and Z-transforms.

.. image:: examples/schematics/lpf1-buffer-loaded2.png
   :width: 11cm

Contents:
=========

.. toctree::
   :maxdepth: 3

   install.rst
   overview.rst
   expressions.rst
   discretetime.rst                         
   circuits.rst
   networks.rst
   netlists.rst
   schematics.rst
   tutorials.rst
   novice.rst              
   latex.rst
   transforms.rst
   config.rst                                                      
   problems.rst                          
   internals.rst
   releases.rst
   modules.rst
   about.rst   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

