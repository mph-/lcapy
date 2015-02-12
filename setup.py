#!/usr/bin/env python

from distutils.core import setup

setup(name = 'lcapy',
    version = '0.6.4-git',
    description = 'Symbolic linear circuit analysis',
    author = 'Michael Hayes',
    requires = [ 'sympy', 'numpy' ],
    author_email = 'michael.hayes@canterbury.ac.nz',
    url = 'https://github.com/mph-/lcapy',
    download_url = 'https://github.com/mph-/lcapy',
    py_modules = [ 'lcapy.core', 'lcapy.netlist', 'lcapy.oneport', 'lcapy.twoport', 'lcapy.threeport', 'lcapy.schematic', 'lcapy.mna', 'lcapy.plot'],
    scripts = ['scripts/schtex.py'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)'
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        ]
    )
