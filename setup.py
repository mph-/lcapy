#!/usr/bin/env python
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='lcapy',
      version='0.39.1',
      author='Michael Hayes',
      author_email='michael.hayes@canterbury.ac.nz',
      description='Symbolic linear circuit analysis',
      long_description = long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/mph-/lcapy',
      download_url='https://github.com/mph-/lcapy',
      # Some of these requirements could be relaxed.
      install_requires=['matplotlib',
                        'scipy',
                        'numpy',                        
                        'sympy',
                        'networkx',
      ],
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
