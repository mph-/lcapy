#!/usr/bin/env python

from distutils.core import setup

setup(name='lcapy',
      version='0.7.1',
      description='Symbolic linear circuit analysis',
      author='Michael Hayes',
      requires=['sympy', 'numpy', 'scipy'],
      author_email='michael.hayes@canterbury.ac.nz',
      url='https://github.com/mph-/lcapy',
      download_url='https://github.com/mph-/lcapy',
      py_modules=['lcapy.core', 'lcapy.netlist', 'lcapy.oneport', 'lcapy.twoport', 'lcapy.threeport', 'lcapy.schematic', 'lcapy.mna', 'lcapy.plot', 'lcapy.latex'],
      scripts=['scripts/schtex.py'],
      license='LGPL'
  )
