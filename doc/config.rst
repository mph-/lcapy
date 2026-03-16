=============
Configuration
=============


Resource parameters
===================

Run-time configuration is controlled by the `rcParams` object.   This is a dictionary of name/value pairs.
The default values can be overriden, for example,

   >>> rcParams['style.schematic'] = 'european'

The default values can also be overridden using an `lcapyrc` file.   This is loaded if it is found in the `.lcapy` directory of the user's home directory.  If an `lcapyrc` file is found in the current working directory it is loaded as well.

Each entry in an `lcapyrc` file is a colon separated name/value pairs.  Entries starting with a `#` are ignored.  Here's an example:

    style.schematic: european


The default entries for `rcParams` are:

    {'sympy.solver': 'DM',
    'sympy.matrix.inverse': 'DM',
    'sympy.print_order': 'lex',
    'functions.heaviside_zero': 0.5,
    'functions.unitstep_zero': 1.0,
    'symbols.imaginary': 'j',
    'symbols.heaviside': 'u',
    'circuit.current_sign_convention': 'passive',
    'schematics.implicit_default': 'sground',
    'schematics.autoground_default': 'sground',
    'schematics.draw_nodes': 'primary',
    'schematics.label_nodes': 'primary',
    'schematics.label_values': True,
    'schematics.label_ids': True,
    'schematics.label_style': 'aligned',
    'schematics.label_value_style': 'eng3',
    'schematics.label_flip': False,
    'schematics.annotate_values': False,
    'schematics.anchor': 'south east',
    'schematics.autoground': 'none',
    'schematics.scale': 1.0,
    'schematics.dpi': 300,
    'schematics.cpt_size': 1.5,
    'schematics.node_spacing': 2.0,
    'schematics.help_lines': 0.0,
    'schematics.style': 'american',
    'schematics.voltage_dir': 'RP',
    'os.linux.ghostscript': 'gs',
    'os.macos.ghostscript': 'gs',
    'os.windows32.ghostscript': 'gswin32c',
    'os.windows64.ghostscript': 'gswin64c',
    'os.linux.convert': 'convert',
    'os.macos.convert': 'convert',
    'os.windows32.convert': 'magick convert',
    'os.windows64.convert': 'magick convert',
    'os.linux.pdf2svg': 'pdf2svg',
    'os.macos.pdf2svg': 'pdf2svg',
    'os.windows32.pdf2svg': 'pdf2svg',
    'os.windows64.pdf2svg': 'pdf2svg',
    'os.linux.pdflatex': 'pdflatex',
    'os.macos.pdflatex': 'pdflatex',
    'os.windows32.pdflatex': 'pdflatex',
    'os.windows64.pdflatex': 'pdflatex',
    'os.linux.dot': 'dot',
    'os.macos.dot': 'dot',
    'os.windows32.dot': 'dot',
    'os.windows64.dot': 'dot',
    'units.loose': True,
    'units.show': False,
    'units.abbreviate': True,
    'units.canonical': False,
    'units.check': True}






Printing
========

There are three ways to print an expression:

1. string (mostly for use in the debugger)
2. unicode (for the python/Ipython command line)
3. latex (for jupyter notebooks and LaTeX files)

For each way, the output representation of an expression can be changed using the dictionaries `str_expr_map`, `pretty_expr_map`, and `latex_expr_map`.  For example, to print the Heaviside "function" as an `H` instead of the default `u`, for the LaTeX representation use:

>>> from lcapy.config import latex_expr_map
>>> latex_expr_map[sym.Heaviside] = 'H'

Handling the the imaginary operator `sym.I` is trickier since this is ingrained as `i` in SymPy.  By default, Lcapy uses `j` but this can be forced to be printed as `i`, see `lcapy/config.py`.


Parsing
=======

Sympy has many predefined symbols.  Lcapy ignores the ones listed in `config.excludes`.  For example, `N` represents the numerical evaluation function.  If you wish to have this symbol, remove it from the `excludes` list, for example,

>>> from lcapy.config import excludes
>>> excludes.remove('N')

Here's another example, that removes the falling factorial function `ff`,

>>> from lcapy.config import excludes
>>> excludes.append('ff')
