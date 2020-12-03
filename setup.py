#!/usr/bin/env python
from setuptools import setup, find_packages


__version__ = '0.72'

with open("README.md", "r") as fh:
    long_description = fh.read()

tests_require = ['nose']


setup(name='lcapy',
      version=__version__,
      author='Michael Hayes',
      author_email='michael.hayes@canterbury.ac.nz',
      description='Symbolic linear circuit analysis',
      long_description = long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/mph-/lcapy',
      download_url='https://github.com/mph-/lcapy',
      install_requires=['matplotlib',
                        'scipy',
                        'numpy',
                        'sympy',
                        'networkx',
                        'setuptools',
                        'wheel'
      ],
      python_requires='>=3.6',
      extras_require={
          'test': tests_require,
          'doc': ['sphinx', 'ipython'],
          'release': ['wheel', 'twine'],
      },
      tests_require=tests_require,
      packages=find_packages(exclude=['demo']),
      entry_points={
          'console_scripts': [
              'schtex=lcapy.scripts.schtex:main',
          ],
      },
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
          "Operating System :: OS Independent",
      ],
)
