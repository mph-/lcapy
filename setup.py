#!/usr/bin/env python

from distutils.core import setup

setup(name='lcapy', version='0.29.1',
      description='Symbolic linear circuit analysis',
      author='Michael Hayes',
      requires=['sympy', 'numpy', 'scipy'],
      author_email='michael.hayes@canterbury.ac.nz',
      url='https://github.com/mph-/lcapy',
      download_url='https://github.com/mph-/lcapy',
      py_modules=['lcapy.core', 'lcapy.netlist', 'lcapy.oneport',
                  'lcapy.twoport', 'lcapy.threeport', 'lcapy.schematic',
                  'lcapy.mna', 'lcapy.plot', 'lcapy.latex', 'lcapy.grammar',
                  'lcapy.parser', 'lcapy.schemcpts', 'lcapy.schemmisc',
                  'lcapy.schemgraph', 'lcapy.mnacpts', 'lcapy.sympify',
                  'lcapy.acdc', 'lcapy.network', 'lcapy.circuit',
                  'lcapy.netfile', 'lcapy.system', 'lcapy.laplace',
                  'lcapy.fourier', 'lcapy.ratfun', 'lcapy.utils', 'lcapy.expr',
                  'lcapy.sexpr', 'lcapy.vector', 'lcapy.matrix',
                  'lcapy.symbols', 'lcapy.cexpr', 'lcapy.texpr',
                  'lcapy.fexpr', 'lcapy.omegaexpr', 'lcapy.sfwexpr',
                  'lcapy.noiseexpr', 'lcapy.phasor', 'lcapy.super',
                  'lcapy.context', 'lcapy.sym', 'lcapy.functions',
                  'lcapy.printing', 'lcapy.config'
      ], scripts=['scripts/schtex.py'],
      license='LGPL' )
