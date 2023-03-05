#!/usr/bin/python3
from os import system
from sys import exit


def main(argv=None):

    system('ipython -c "from lcapy import *" -i')


if __name__ == '__main__':
    exit(main())
