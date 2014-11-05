#!/usr/bin/python
"""schtex V0.1
Copyright (c) 2014 Michael P. Hayes, UC ECE, NZ

Usage: schtex file.sch
"""

from optparse import OptionParser
import sys
from lcapy import Schematic


def main (argv=None):

    if argv is None:
        argv = sys.argv

    version = __doc__.split('\n')[0]

    parser = OptionParser(usage = '%prog schematic-file', version = version, 
                          description = __doc__)

    parser.add_option('--label-nodes', action = 'store_true',
                      dest = 'label_nodes', default = True,
                      help = 'label nodes')
    
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error('missing argument')
        return 1

    sch = Schematic(args[0])

    sch.draw(label_nodes=options.label_nodes)

    return 0


if __name__ == '__main__':
    sys.exit(main())
