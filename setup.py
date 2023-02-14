#!/usr/bin/env python
from setuptools import setup, find_packages


__version__ = '1.12'

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

tests_require = ['nose', 'flake8', 'flake8-bugbear',
                 'flake8-comprehensions', 'flake8-requirements']


setup(name='lcapy',
      version=__version__,
      author='Michael Hayes',
      author_email='michael.hayes@canterbury.ac.nz',
      description='Symbolic linear circuit analysis',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/mph-/lcapy',
      download_url='https://github.com/mph-/lcapy',
      install_requires=['matplotlib',
                        'scipy',
                        'numpy',
                        'sympy>=1.10.1',  # >=1.7.1 works
                        'networkx',
                        'IPython',
                        'setuptools',
                        'wheel',
                        'property_cached'
                        ],
      python_requires='>=3.7',  # >=3.6 should still work but not tested
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
              'lcapy=lcapy.scripts.lcapy:main',
          ],
      },
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
          "Operating System :: OS Independent",
      ],
      )
