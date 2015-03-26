#!/usr/bin/python
"""schtex V0.1
Copyright (c) 2014 Michael P. Hayes, UC ECE, NZ

Usage: schtex infile.sch [outfile.tex|pdf|png|svg]
"""

from __future__ import print_function
from optparse import OptionParser
import sys
import os


def main (argv=None):

    if argv is None:
        argv = sys.argv

    version = __doc__.split('\n')[0]

    parser = OptionParser(usage='%prog schematic-file [output-file]', version=version, 
                          description=__doc__)

    parser.add_option('--label-nodes', type='str',
                      dest='label_nodes', default=None,
                      help='label nodes, choice none, primary, all')

    parser.add_option('--nolabel-nodes', action='store_false',
                      dest='label_nodes',
                      help="don't label nodes")

    parser.add_option('--draw-nodes', type='str',
                      dest='draw_nodes', default=None,
                      help='draw nodes, choice none, primary, all')

    parser.add_option('--nodraw-nodes', action='store_false',
                      dest='draw_nodes',
                      help="don't draw nodes")

    parser.add_option('--s-model', action='store_true',
                      dest='s_model', default=False,
                      help='generate s-domain model schematic')

    parser.add_option('--p-model', action='store_true',
                      dest='p_model', default=False,
                      help='generate pre-initial model schematic')

    parser.add_option('--k-model', action='store_true',
                      dest='k_model', default=False,
                      help='generate schematic with independent sources killed')

    parser.add_option('--scale', type='float',
                      dest='scale', default=1,
                      help='schematic scale factor')

    parser.add_option('--stretch', type='float',
                      dest='stretch', default=1,
                      help='schematic stretch factor')
    
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error('missing argument')
        return 1

    infilename = args[0]
    outfilename = None
    if len(args) > 1:
        outfilename = args[1]

    from lcapy import Circuit

    cct = Circuit(infilename)
    if options.k_model:
        cct = cct.kill()
    if options.s_model:
        cct = cct.s_model()
    if options.p_model:
        cct = cct.pre_initial_model()

    if options.label_nodes not in ('none', 'all', 'primary', False, None):
        raise ValueError('Illegal option %s for label_nodes' % options.label_nodes)

    if options.draw_nodes not in ('none', 'all', 'primary', False, None):
        raise ValueError('Illegal option %s for draw_nodes' % options.draw_nodes)

    cct.draw(label_nodes=options.label_nodes, draw_nodes=options.draw_nodes,
             filename=outfilename, scale=options.scale, stretch=options.stretch)

    return 0


if __name__ == '__main__':
    sys.exit(main())
