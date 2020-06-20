.. _installation:

============
Installation
============


Requirements
============

Lcapy requires the following python packages: python, lcapy, sympy, numpy, matplotlib, networkx

For schematic drawing you require:

1. pdflatex

2. circuitikz (https://www.ctan.org/pkg/circuitikz).  Lcapy currently
   expects circuitikz version 1.0.1.  Unfortunately, different
   circuitikz releases tweak component sizes.

3. ghostscript (http://www.ghostscript.com/), imagemagick convert (http://www.imagemagick.org/), or `pdftoppm`

For nice maths formatting in a jupyter notebook you require mathjax (https://www.mathjax.org/).  This is not essential; if it is not loaded then an active internet connection is required.


Installation for Linux
======================

1. Lcapy (and its python dependencies) can be installed using:

.. code-block:: console
                
    $ pip3 install lcapy

2. For systems using apt, all the other packages required for Lcapy can be installed using:

.. code-block:: console
                
    $ sudo apt install texlive-latex-base texlive-pictures texlive-latex-extra imagemagick ghostscript libjs-mathjax fonts-mathjax

3. Download circuitikz (https://www.ctan.org/pkg/circuitikz)


Installation for Windows
========================

1. Lcapy (and its python dependencies) can be installed using `pip3 install lcapy`

2. To install ghostscript see https://www.ghostscript.com/download/gsdnld.html
   
3. Download circuitikz (https://www.ctan.org/pkg/circuitikz)
      

Installation for development
============================

1. Install Lcapy and its dependencies as per instructions for your platform

2. Download Lcapy sources from https://github.com/mph-/lcapy as a .zip file or preferably using git:

.. code-block:: console
                     
   $ git clone https://github.com/mph-/lcapy

3.  Install Lcapy using:

.. code-block:: console
    
   $ cd lcapy
   $ sudo python setup.py install

Note, if you do not have root access, you can use  `virtualenv` or  set the environment variable `PYTHONPATH` to find the source files for Lcapy.

4. For testing you need `nosetests`.  For example, using apt:

.. code-block:: console
                     
   $ sudo apt install python3-nose


5. For building the docs you need `sphinx`.  For example, using apt:

.. code-block:: console
                     
   $ sudo apt install python3-sphinx
   
6. For making releases you need `wheel` and `twine`

   $ sudo apt install python3-wheel
   $ pip3 install twine   

7. For debugging schematic graphs `dot` is required:

.. code-block:: console
                     
   $ sudo apt install graphviz   

  
