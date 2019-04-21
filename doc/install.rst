.. _installation:

============
Installation
============

Lcapy (and its python dependencies) can be installed using:

>>> pip3 install lcapy

For schematic drawing you require:

1. pdflatex

2. circuitikz (https://www.ctan.org/pkg/circuitikz).  Lcapy currently
   expects circuitikz version 2017/05/28.  Unfortunately, different
   circuitikz releases tweak component sizes.

3. imagemagick convert (http://www.imagemagick.org/) or `pdftoppm`

4. ghostscript (http://www.ghostscript.com/)

For nice maths formatting in a jupyter notebook you require mathjax (https://www.mathjax.org/).  This is not essential; if it is not loaded then an active internet connection is required.


Installation for Linux (Ubuntu and variants)
============================================

All the other packages required for Lcapy can be installed using:

   >>> sudo apt install texlive-latex-base texlive-pictures texlive-latex-extra imagemagick ghostscript libjs-mathjax fonts-mathjax

   
Installation from github
========================

1. Lcapy can be downloaded from https://github.com/mph-/lcapy as a .zip file or preferably using git::
     
   >>> git clone https://github.com/mph-/lcapy

2.  You will also need to install scipy, numpy, matplotlib, and networkx.
   
3.  Lcapy can be installed using:

  >>> cd lcapy
  >>> sudo python setup.py install

  Note, if you do not have root access, you can use  `virtualenv` or  set the environment variable `PYTHONPATH` to find the source files for Lcapy.
