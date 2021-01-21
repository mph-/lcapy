.. _installation:

============
Installation
============


Requirements
============

Lcapy requires the following Python packages: scipy, sympy, numpy, matplotlib, networkx

For schematic drawing you require:

1. pdflatex

2. CircuiTikZ (https://www.ctan.org/pkg/circuitikz).  Lcapy currently
   expects CircuiTikZ version 1.0.1.  Unfortunately, different
   CircuiTikZ releases tweak component sizes.

3. Ghostscript (http://www.ghostscript.com/), ImageMagick (http://www.imagemagick.org/), or `pdftoppm`

For nice maths formatting in a Jupyter notebook you require mathjax (https://www.mathjax.org/).  This is not essential; if it is not loaded then an active internet connection is required.


Installation for Linux
======================

1. Lcapy (and its Python dependencies) can be installed using:

.. code-block:: console

    $ pip3 install lcapy

2. For systems using apt, all the other packages required for Lcapy can be installed using:

.. code-block:: console

    $ sudo apt install ghostscript texlive-latex-extra


Installation for macOS
======================

1. Lcapy (and its Python dependencies) can be installed using:

.. code-block:: console

    $ pip3 install lcapy

2. Using Homebrew (https://brew.sh), all the other packages required for Lcapy can be installed using:

.. code-block:: console

    $ brew install mactex

or to minimise disk space, install a smaller set of LaTeX packages using:

.. code-block:: console

    $ brew install basictex
    $ eval "$(/usr/libexec/path_helper)"
    $ sudo tlmgr install collection-latexextra


Installation for Windows
========================

1. Lcapy (and its Python dependencies) can be installed using

.. code-block:: console

   $ pip3 install lcapy`

2. To install ghostscript see https://www.ghostscript.com/download/gsdnld.html

3. Download circuitikz (https://www.ctan.org/pkg/circuitikz)


Installation for development
============================

1. Install Lcapy's dependencies as per instructions for your platform

2. Download Lcapy sources from https://github.com/mph-/lcapy as a .zip file or preferably using git:

.. code-block:: console

   $ git clone https://github.com/mph-/lcapy

3.  Install Lcapy using:

.. code-block:: console
                
   $ cd lcapy
   $ pip3 install --editable .[test,release]

4. For building the docs you need `sphinx`, `ipython` and `pycairo`.  For example, using apt:

.. code-block:: console

   $ sudo apt-get install gir1.2-gtk-3.0 python3-gi python3-gi-cairo
   $ pip3 install --editable .[doc]

5. For debugging schematic graphs `dot` is required:

.. code-block:: console

   $ sudo apt install graphviz

6. To run style guide checking locally

.. code-block:: console

   $ pip3 install flake8 flake8-bugbear flake8-requirements flake8-comprehensions
   
7. For coverage analysis

.. code-block:: console

   $ pip3 install coverage
   
