#!/usr/bin/env python
"""schtex V0.1
Copyright (c) 2014--2019 Michael P. Hayes, UC ECE, NZ

Usage: schtex infile.sch [outfile.tex|pdf|png|svg|sch|pgf]
"""

from __future__ import print_function
from optparse import OptionParser
import sys
import os

def schtex_exception(type, value, tb):
   if hasattr(sys, 'ps1') or not sys.stderr.isatty():
      # We are not in interactive mode or we don't have a tty-like
      # device, so call the default hook
      sys.__excepthook__(type, value, tb)
   else:
      import traceback, pdb
      # We are in interactive mode, print the exception...
      traceback.print_exception(type, value, tb)
      print()
      # ...then start the debugger in post-mortem mode.
      pdb.pm()

def main (argv=None):

    if argv is None:
        argv = sys.argv

    version = __doc__.split('\n')[0]

    parser = OptionParser(usage="""
%prog schematic-file [outfile.tex|pdf|png|svg|sch|pgf]""",
                          version=version, description=__doc__)

    parser.add_option('--debug', action='store_true',
                      dest='debug', default=None,
                      help="enable debugging")

    parser.add_option('--draw-nodes', type='str',
                      dest='draw_nodes', default=None,
                      help='draw nodes, choice: none, connections, primary, all')

    parser.add_option('--nodraw-nodes', action='store_false',
                      dest='draw_nodes',
                      help="don't draw nodes")

    parser.add_option('--label-nodes', type='str',
                      dest='label_nodes', default=None,
                      help='label nodes, choice: none, alpha, pins, primary, all or a list of comma separated node names in braces, e.g, {1, 2, 3}')

    parser.add_option('--nolabel-nodes', action='store_false',
                      dest='label_nodes',
                      help="don't label nodes")

    parser.add_option('--nolabel-ids', action='store_false',
                      dest='label_ids', default=None,
                      help="don't label element identifiers")

    parser.add_option('--label-ids', action='store_true',
                      dest='label_ids', default=None,
                      help="label element identifiers")

    parser.add_option('--nolabel-values', action='store_false',
                      dest='label_values', default=None,
                      help="don't label values")

    parser.add_option('--label-values', action='store_true',
                      dest='label_values', default=None,
                      help="labels values")

    parser.add_option('--s-model', action='store_true',
                      dest='s_model', default=False,
                      help='generate s-domain model schematic')

    parser.add_option('--ac-model', action='store_true',
                      dest='ac_model', default=False,
                      help='generate AC model schematic')

    parser.add_option('--p-model', action='store_true',
                      dest='p_model', default=False,
                      help='generate pre-initial model schematic')

    parser.add_option('--k-model', action='store_true',
                      dest='k_model', default=False,
                      help='generate schematic with independent sources killed')

    parser.add_option('--scale', type='float',
                      dest='scale', default=None,
                      help='schematic scale factor, this scales the schematic size but not the fonts')

    parser.add_option('--node-spacing', type='float',
                      dest='node_spacing', default=None,
                      help='this specifies the spacing of the nodes of a component')

    parser.add_option('--cpt-size', type='float',
                      dest='cpt_size', default=None,
                      help='this specifies the size of a component; it needs to be smaller than node_spacing')
    
    parser.add_option('--help-lines', type='float',
                      dest='help_lines', default=None,
                      help="draw help lines")

    parser.add_option('--xgraph', action='store_true',
                      dest='xgraph', default=False,
                      help="generate graph of component x positions")

    parser.add_option('--ygraph', action='store_true',
                      dest='ygraph', default=False,
                      help="generate graph of component y positions")

    parser.add_option('--stage', type='int',
                      dest='stage', default=0,
                      help='graph analysis stage')

    parser.add_option('--pdb', action='store_true',
                      default=False,
                      help="enter python debugger on exception")    

    parser.add_option('--renumber', type='str',
                      dest='renumber', default=None,
                      help='renumber nodes, e.g, 10:1, 11:2')
    
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error('missing argument')
        return 1

    infilename = args[0]
    outfilename = None
    if len(args) > 1:
        outfilename = args[1]

    if options.pdb:
        sys.excepthook = schtex_exception
        
    from lcapy import Circuit

    cct = Circuit(infilename)
    if options.k_model:
        cct = cct.kill()
    if options.s_model:
        cct = cct.s_model()
    if options.ac_model:
        cct = cct.ac_model()
    if options.p_model:
        cct = cct.pre_initial_model()

    if options.renumber:
        parts = options.renumber.split(',')
        node_map = {}
        for part in parts:
            part = part.strip()
            fields = part.split(':')
            if len(fields) != 2:
                raise ValueError('Expecting mapping of form a:b got %s' % part)
            node_map[fields[0]] = fields[1]
        cct = cct.renumber(node_map)

    if options.label_nodes not in ('none', 'all', 'alpha', 'pins', 'primary', False, None):
        raise ValueError('Illegal option %s for label_nodes' % options.label_nodes)

    if options.draw_nodes not in ('none', 'all', 'primary', 'connections',
                                  False, None):
        raise ValueError('Illegal option %s for draw_nodes' % options.draw_nodes)

    if outfilename is not None:
       base, ext = os.path.splitext(outfilename)
       if ext == '.sch':
          open(outfilename, 'w').write(str(cct))
          return 0
        
    nosave = options.xgraph or options.ygraph

    if not options.xgraph and not options.ygraph:
        cct.draw(label_nodes=options.label_nodes,
                 draw_nodes=options.draw_nodes,
                 label_ids=options.label_ids,
                 label_values=options.label_values, 
                 filename=outfilename, scale=options.scale,
                 node_spacing=options.node_spacing, cpt_size=options.cpt_size,
                 help_lines=options.help_lines, debug=options.debug,
                 nosave=nosave)

    if options.xgraph:
        cct.sch.make_graphs()
        cct.sch.xgraph.dot(outfilename, stage=options.stage)
        return 0

    if options.ygraph:
        cct.sch.make_graphs()
        cct.sch.ygraph.dot(outfilename, stage=options.stage)
        return 0

    return 0


if __name__ == '__main__':
    sys.exit(main())
