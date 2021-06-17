
.. lcapy documentation master file, created by
   sphinx-quickstart on Mon Feb 17 17:08:01 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Lcapy's documentation!
=================================

Lcapy is a Python package for symbolic linear circuit analysis and signal-processing.  It uses SymPy for the underlying symbolic analysis and infrastructure.  In addition, Lcapy can semi-automate the drawing of high-quality schematics from a netlist, including diodes, transistors, and other non-linear components.  It can also symbolically transform both continuous-time and discrete-time signals, using Laplace, Fourier, and z-transforms.  Expressions can be formatted in many forms and plotted with auto-labelling of the axes.

.. image:: examples/schematics/lpf1-buffer-loaded2.png
   :width: 11cm

_ 
   >>> H = 5 * (s**2 + 1) / (s**2 + 5*s + 4)
   >>> H(t)   
                -t       -4⋅t           
            10⋅ℯ     85⋅ℯ               
   5⋅δ(t) + ────── - ────────  for t ≥ 0
              3         3               

   
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

