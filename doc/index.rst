
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

If you have suggestions for improvements, please create an issue at
https://github.com/mph-/lcapy/issues.

For a graphical user interface, see https://github.com/mph-/lcapy-gui.

If you are using Lcapy for teaching, please email to `''.join(reversed(['nz', '.', 'ac', '.', 'canterbury', '@', 'hayes', '.', 'michael']))`.

To cite Lcapy in publications use `Hayes M. 2022. Lcapy: symbolic linear circuit analysis with Python. PeerJ Computer Science 8:e875` https://doi.org/10.7717/peerj-cs.875


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
   systems.rst
   schematics.rst
   gallery.rst
   tutorials.rst
   novice.rst
   advanced.rst
   latex.rst
   transforms.rst
   config.rst
   applications.rst
   problems.rst
   internals.rst
   releases.rst
   development.rst
   modules.rst
   about.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
