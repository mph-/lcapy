.. _installation:

============
Installation
============

- You will need to install SymPy, see http://docs.sympy.org/latest/install.html.

- For plotting you will need to install matplotlib (this requires Numpy).

- Lcapy can be downloaded from https://github.com/mph-/lcapy using git or as a .zip file.


Lcapy installation for Linux
============================

- You need to install the Sympy package:

  >>> sudo apt-get install python-sympy

- For plotting you need the Matplotlib and Numpy packages:

  >>> sudo apt-get install python-numpy python-matplotlib

- The easiest way to obtain the Lcapy sources is to use git:

  >>> git clone https://github.com/mph-/lcapy

- Lcapy can be installed using:

  >>> cd lcapy
  >>> sudo python setup.py install

- If you do not have root access, you can set the environment variable `PYTHONPATH` to find the source files for Lcapy.

- While you are at it, it is worthwhile to install ipython, the interactive python shell

  >>> sudo apt-get install ipython
