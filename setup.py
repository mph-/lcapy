#!/usr/bin/env python

from distutils.core import setup

setup(name = 'lcapy',
    version = '0.2.0',
    description = 'Linear circuit analysis',
    author = 'Michael Hayes',
    requires = [ 'sympy', 'numpy' ],
    author_email = 'michael.hayes@canterbury.ac.nz',
    url = 'https://github.com/mph-/lcapy',
    download_url = 'https://github.com/mph-/lcapy',
    py_modules = [ 'lcapy.mcircuit', 'lcapy.netlist' ],
    scripts = [ ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)'
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        ]
    )
