#!/usr/bin/env python
"""schtex V0.2
Copyright (c) 2014--2020 Michael P. Hayes, UC ECE, NZ

Usage: schtex infile.sch [outfile.tex|pdf|png|svg|sch|pgf]
"""

from __future__ import print_function
from matplotlib.pyplot import show
from argparse import ArgumentParser
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

    parser = ArgumentParser(description='Convert Lcapy schematic to image file.')
    parser.add_argument('--version', action='version', version=__doc__.split('\n')[0])
    
    parser.add_argument('--debug', action='store_true',
                        dest='debug', default=None,
                        help="enable debugging")
    
    parser.add_argument('--draw-nodes', type=str,
                        dest='draw_nodes', default=None,
                        help='draw nodes, choice: none, connections, primary, all')
    
    parser.add_argument('--nodraw-nodes', action='store_false',
                        dest='draw_nodes',
                        help="don't draw nodes")
    
    parser.add_argument('--label-nodes', type=str,
                        dest='label_nodes', default=None,
                        help='label nodes, choice: none, alpha, pins, primary, all or a list of comma separated node names in braces, e.g, {1, 2, 3}')
    
    parser.add_argument('--nolabel-nodes', action='store_false',
                        dest='label_nodes',
                        help="don't label nodes")
    
    parser.add_argument('--nolabel-ids', action='store_false',
                        dest='label_ids', default=None,
                        help="don't label element identifiers")
    
    parser.add_argument('--label-ids', action='store_true',
                        dest='label_ids', default=None,
                        help="label element identifiers")
    
    parser.add_argument('--nolabel-values', action='store_false',
                        dest='label_values', default=None,
                        help="don't label values")
    
    parser.add_argument('--label-values', action='store_true',
                        dest='label_values', default=None,
                        help="labels values")
    
    parser.add_argument('--s-model', action='store_true',
                        dest='s_model', default=False,
                        help='generate s-domain model schematic')

    parser.add_argument('--noisy-model', action='store_true',
                        dest='noisy_model', default=False,
                        help='generate noisy model schematic')
    
    parser.add_argument('--ac-model', action='store_true',
                        dest='ac_model', default=False,
                        help='generate AC model schematic')
    
    parser.add_argument('--p-model', action='store_true',
                        dest='p_model', default=False,
                        help='generate pre-initial model schematic')
    
    parser.add_argument('--k-model', action='store_true',
                        dest='k_model', default=False,
                        help='generate schematic with independent sources killed')

    parser.add_argument('--kill-zero', action='store_true',
                        default=False,
                        help='generate schematic with zero independent sources killed')    
    
    parser.add_argument('--scale', type=float,
                        dest='scale', default=None,
                        help='schematic scale factor, this scales the schematic size but not the fonts')
    
    parser.add_argument('--node-spacing', type=float,
                        dest='node_spacing', default=None,
                        help='this specifies the spacing of the nodes of a component')
    
    parser.add_argument('--cpt-size', type=float,
                        dest='cpt_size', default=None,
                        help='this specifies the size of a component; it needs to be smaller than node_spacing')

    parser.add_argument('--dpi', type=float, default=None,
                        help='this specifies the number of dots per inch when creating a png file')    
    
    parser.add_argument('--help-lines', type=float,
                        dest='help_lines', default=None,
                        help="draw help lines")
    
    parser.add_argument('--xgraph', action='store_true',
                        dest='xgraph', default=False,
                        help="generate graph of component x positions")

    parser.add_argument('--circuitgraph', action='store_true',
                        default=False,
                        help="generate circuit graph")    
    
    parser.add_argument('--ygraph', action='store_true',
                        dest='ygraph', default=False,
                        help="generate graph of component y positions")
    
    parser.add_argument('--stage', type=int,
                        dest='stage', default=0,
                        help='select graph analysis stage 0--4')
    
    parser.add_argument('--pdb', action='store_true',
                        default=False,
                        help="enter python debugger on exception")

    parser.add_argument('--allow-anon', action='store_true',
                        default=False,
                        help="allow anonymous component ids")        
    
    parser.add_argument('--renumber', type=str,
                        dest='renumber', default=None,
                        help='renumber nodes, e.g, 10:1, 11:2 or all')

    parser.add_argument('--font', type=str, default=None,
                        help=r'specify default font, e.g., \small, \sffamily\tiny')    

    parser.add_argument('--options', type=str, default=None,
                        help='specify circuitikz options')
    
    parser.add_argument('--preamble', type=str, default=None,
                        help='specify circuitikz commands at start')

    parser.add_argument('--postamble', type=str, default=None,
                        help='specify circuitikz commands at end')

    parser.add_argument('filename', type=str, help='schematic filename')
    parser.add_argument('outfilename', type=str, nargs='?',
                        help='output filename', default=None)    

    parser.add_argument('--show', action='store_true',
                        dest='show', default=None,
                        help="show image")
    
    args = parser.parse_args()    

    infilename = args.filename
    outfilename = args.outfilename

    if args.pdb:
        sys.excepthook = schtex_exception
        
    from lcapy import Circuit, CircuitGraph

    cct = Circuit(infilename, allow_anon=args.allow_anon)
    if args.k_model:
        cct = cct.kill()
    if args.kill_zero:
        cct = cct.kill_zero()        
    if args.s_model:
        cct = cct.s_model()
    if args.ac_model:
        cct = cct.ac_model()
    if args.p_model:
        cct = cct.pre_initial_model()
    if args.noisy_model:
        cct = cct.noisy()        

    if args.renumber == 'all':
       cct = cct.renumber()
    elif args.renumber:
        parts = args.renumber.split(',')
        node_map = {}
        for part in parts:
            part = part.strip()
            fields = part.split(':')
            if len(fields) != 2:
                raise ValueError('Expecting mapping of form a:b got %s' % part)
            node_map[fields[0]] = fields[1]
        cct = cct.renumber(node_map)

    if args.label_nodes not in ('none', 'all', 'alpha', 'pins', 'primary', False, None):
        raise ValueError('Illegal option %s for label_nodes' % args.label_nodes)

    if args.draw_nodes not in ('none', 'all', 'primary', 'connections',
                                  False, None):
        raise ValueError('Illegal option %s for draw_nodes' % args.draw_nodes)

    if outfilename is not None:
       base, ext = os.path.splitext(outfilename)
       if ext == '.sch':
          open(outfilename, 'w').write(str(cct))
          return 0
        
    nosave = args.xgraph or args.ygraph

    if not args.xgraph and not args.ygraph and not args.circuitgraph:
        kwargs = {}
        if args.options != None:
            kwargs['options'] = args.options
        if args.preamble != None:
            kwargs['preamble'] = args.preamble
        if args.postamble != None:
            kwargs['postamble'] = args.postamble            
        if args.font != None:
            kwargs['font'] = args.font          
       
        cct.draw(label_nodes=args.label_nodes,
                 draw_nodes=args.draw_nodes,
                 label_ids=args.label_ids,
                 label_values=args.label_values, 
                 filename=outfilename, scale=args.scale,
                 node_spacing=args.node_spacing, cpt_size=args.cpt_size,
                 help_lines=args.help_lines, debug=args.debug,
                 dpi=args.dpi, nosave=nosave, **kwargs)

    if args.xgraph:
        cct.sch.make_graphs(debug=args.debug)
        cct.sch.xgraph.dot(outfilename, stage=args.stage)

    if args.ygraph:
        cct.sch.make_graphs(debug=args.debug)
        cct.sch.ygraph.dot(outfilename, stage=args.stage)

    if args.circuitgraph:
       cg = CircuitGraph(cct)
       cg.draw(filename=outfilename)

    if args.show:
        show()
        
    return 0


if __name__ == '__main__':
    sys.exit(main())
