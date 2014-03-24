#!/usr/bin/env python

from distutils.core import setup

setup(name = 'mcircuit',
    version = '0.1.0',
    description = 'Linear circuit analysis',
    author = 'Michael Hayes',
    requires = [ 'sympy', 'numpy' ],
    author_email = 'michael.hayes@canterbury.ac.nz',
    url = 'https://github.com/mph-/mcircuit',
    download_url = 'https://github.com/mph-/mcircuit',
    py_modules = [ 'mcircuit' ],
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
