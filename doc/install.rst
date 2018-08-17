.. _installation:

============
Installation
============

- You need to install SymPy, see http://docs.sympy.org/latest/install.html.

- For plotting you need to install matplotlib (this requires NumPy).

- For schematic drawing you need pdflatex, the circuitikz package (https://www.ctan.org/pkg/circuitikz), the imagemagick convert program (http://www.imagemagick.org/), and ghostscript (http://www.ghostscript.com/).

- For nice rendering of maths you need mathjax
  (https://www.mathjax.org/).  This is not essential, if it is not
  loaded then an active internet connection is required.

- Lcapy can be downloaded from https://github.com/mph-/lcapy using git or as a .zip file.


Lcapy installation for Linux (Ubuntu and variants)
==================================================

The following instructions are for Python2.7.  For Python 3, replace
the python with python3 in the package names and when running setup.py.

- You need to install the SymPy package:

  >>> sudo apt install python-sympy

- For plotting you need the Matplotlib and NumPy packages:

  >>> sudo apt install python-numpy python-matplotlib

- For schematic drawing you need pdflatex, the circuitikz package,
  imagemagick, and ghostscript.

  >>> sudo apt install texlive-latex-base texlive-pictures texlive-latex-extra imagemagick ghostscript

- For nice rendering of maths you need mathjax.  This is not
  essential, if it is not loaded then an active internet connection is
  required.

  >>> sudo apt install libjs-mathjax fonts-mathjax

- The easiest way to obtain the Lcapy sources is to use git:

  >>> git clone https://github.com/mph-/lcapy

- Lcapy can be installed using:

  >>> cd lcapy
  >>> sudo python setup.py install

- If you do not have root access, you can set the environment variable `PYTHONPATH` to find the source files for Lcapy.

- While you are at it, it is worthwhile to install ipython, the interactive python shell (this is also useful for displaying notebooks in a web browser).

  >>> sudo apt install ipython

