.. _installation:

============
Installation
============

- You need to install SymPy, see http://docs.sympy.org/latest/install.html.

- For plotting you need to install matplotlib (this requires Numpy).

- For schematic drawing you need pdflatex and the circuitikz package.

- Lcapy can be downloaded from https://github.com/mph-/lcapy using git or as a .zip file.


Lcapy installation for Linux
============================

- You need to install the Sympy package:

  >>> sudo apt-get install python-sympy

- For plotting you need the Matplotlib and Numpy packages:

  >>> sudo apt-get install python-numpy python-matplotlib

- For schematic drawing you need pdflatex and the circuitikz package.

- The easiest way to obtain the Lcapy sources is to use git:

  >>> git clone https://github.com/mph-/lcapy

- Lcapy can be installed using:

  >>> cd lcapy
  >>> sudo python setup.py install

- If you do not have root access, you can set the environment variable `PYTHONPATH` to find the source files for Lcapy.

- While you are at it, it is worthwhile to install ipython, the interactive python shell  (this is also useful for displaying notebooks in a web browser)

  >>> sudo apt-get install ipython
